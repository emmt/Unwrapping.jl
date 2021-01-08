module Unwrapping

using FFTW

"""
    twopi(T)

yields the value of `2π` with floating-point type `T`.

"""
twopi(T::Type{<:AbstractFloat}) = 2*T(pi)

"""
    arc(x)

yields `x` wrapped in the range `[-π,+π]`.

"""
arc(x::Real) = arc(Float64(x))
arc(x::T) where {T<:AbstractFloat} = wrap(x, twopi(T))
# NOTE: Using `wrap(x, 2π)` is faster than `mod2pi(x + π) - π`.

"""
    wrap(x, p)

yields `x` wrapped in the range `[-p/2,+p/2]`, argument `p` is the period.  The
result is of floating-point type.

"""
wrap(x::Real, p::Real) = wrap(promote(x, p)...)
wrap(x::T, p::T) where {T<:Real} = wrap(Float64(x), Float64(p))
wrap(x::T, p::T) where {T<:AbstractFloat} = x - p*round(x/p)

"""
    unwrap(psi, p=2π; debiased=true) -> phi

yields the unwrapped phase `phi` that best fits, in a least squares sense, the
wrapped phase `psi`.  Optional argument `p` is the period of the wrapping.
The result has an undetermined constant offset.

The method is based on algorithm *"Robust two-dimensional weighted and
unweighted phase unwrapping that uses fast transforms and iterative methods"*
by Dennis C. Ghiglia and Louis A. Romero (J. Opt. Soc. Am. A, Vol. 11,
pp. 107-117, 1994).

Keyword `debiased` (true by default) specifies whether a constant bias (modulo
the period) should be avoided.  Original algorithm yields a biased solution up
to an undetermined additive constant.  Estimating this bias has a cost, so use
`debiased = false` if such an error is irrelevant for your needs.  In that
case, a zero-mean solution is returned (up to rounding errors).

"""
function unwrap(psi::AbstractMatrix{T},
                period::Real; kwds...) where {T<:AbstractFloat}
    unwrap(psi, T(period); kwds...)
end

function unwrap(psi::AbstractMatrix{T},
                period::T = twopi(T);
                debiased::Bool = true) where {T<:AbstractFloat}
    Base.has_offset_axes(psi) && error("input array has non-standard axes")
    n1, n2 = size(psi)

    # We need 2 workspace arrays.
    wrk1 = Array{T}(undef, n1, n2)
    wrk2 = Array{T}(undef, n1, n2)

    # Compute wrapped phase Laplacian by finite differences.  We first compute
    # wrapped phase gradient along a given dimension by finite differences and
    # according to Eq. (2) or Eq. (3).  We then compute the curvature along
    # each dimension by finite differences of the gradient (in-place).  We
    # finally adds the curvatures to form the Laplacian.
    d1 = wrk1
    @inbounds for j2 in 1:n2
        @simd for j1 in 1:n1-1
            d1[j1,j2] = wrap(psi[j1+1,j2] - psi[j1,j2], period)
        end
        let j1 = n1
            d1[j1,j2] = 0
        end
        @simd for j1 in n1:-1:2
            d1[j1,j2] -= d1[j1-1,j2]
        end
    end
    d2 = wrk2
    @inbounds for j2 in 1:n2-1
        @simd for j1 in 1:n1
            d2[j1,j2] = wrap(psi[j1,j2+1] - psi[j1,j2], period)
        end
    end
    let j2 = n2
        @inbounds @simd for j1 in 1:n1
            d2[j1,j2] = 0
        end
    end
    @inbounds for j2 in n2:-1:2
        @simd for j1 in 1:n1
            d2[j1,j2] -= d2[j1,j2-1]
        end
    end
    rho = wrk1 # we re-use one of the workspaces to spare allocations
    @inbounds @simd for j1 in eachindex(wrk1, wrk2)
        rho[j1] = d1[j1] + d2[j1]
    end

    # Compute Q, the inverse of the DCT filter, that is the denominator in
    # Eq. (13) but set Q to 1 at the zero-th frequency to avoid division by
    # zero.  The zero-th frequency being undetermined is a direct consequence
    # of the phase being reconstructed from phase diffeences, it known up to a
    # constant offset.
    q = wrk2 # we re-use the other workspace to spare allocations
    q1 = Array{T}(undef, n1) # array of cosine along first dimension
    a1 = T(π)/n1
    @inbounds @simd for j1 in 1:n1
        q1[j1] = two_cos_minus_two(a1*(j1 - 1))
    end
    if n1 == n2
        @inbounds for j2 in 1:n2
            q2 = q1[j2]
            @simd for j1 in 1:n1
                q[j1,j2] = q1[j1] + q2
            end
        end
    else
        a2 = T(π)/n2
        @inbounds for j2 in 1:n2
            q2 = two_cos_minus_two(a2*(j2 - 1))
            @simd for j1 in 1:n1
                q[j1,j2] = q1[j1] + q2
            end
        end
    end
    q[1,1] = 1 # to avoid division by zero

    # Compute the solution (Julia's dct is of type II).
    dct!(rho)
    rho[1,1] = 0 # should be exactly zero without rounding errors
    phi = rho # we re-use the same workspace to spare allocations
    @inbounds @simd for j in eachindex(phi, rho, q)
        phi[j] = rho[j]/q[j]
    end
    idct!(phi)

    # Fix the piston (modulo the period) in the computed solution.
    debiased && fix_bias!(phi, psi, period)

    return phi
end

"""
    fft_unwrap(psi, p=2π; debiased=true) -> phi

yields the same result as `unwrap` but using a Fast Fourier Transform (FFT)
instead of a Cosine Fourier Transform (DCT).  This version is slower then the
one using the DCT and is provided as an example of how to use the FFT instead
of the DCT.

"""
function fft_unwrap(psi::AbstractMatrix{T},
                    period::T = twopi(T);
                    debiased::Bool = true) where {T<:AbstractFloat}
    Base.has_offset_axes(psi) && error("input array has non-standard axes")
    n1, n2 = size(psi)

    # We need 2 workspace arrays.
    wrk1 = Array{T}(undef, n1, n2)
    wrk2 = Array{T}(undef, n1, n2)

    # The mean level and mean slope (piston and tip-tilt) are undetermined when
    # the FFT is used instead of the DCT.  We therefore compute the discrete
    # wrapped gradient with the mean slope subtracted and then compute the
    # discrete curvature as before to estimate the discrete Laplacian.  After
    # reconstruction, we add the mean slope and the estimated mean level
    # (modulo the period).

    # Fist pass to estimate the mean wrapped slope along each dimension.
    s1 = zero(T)
    if n1 ≥ 2 && n2 ≥ 1
        @inbounds for j2 in 1:n2
            @simd for j1 in 1:n1-1
                s1 += wrap(psi[j1+1,j2] - psi[j1,j2], period)
            end
        end
        s1 /= (n1 - 1)*n2
    end
    s2 = zero(T)
    if n1 ≥ 1 && n2 ≥ 2
        @inbounds for j2 in 1:n2-1
            @simd for j1 in 1:n1
                s2 += wrap(psi[j1,j2+1] - psi[j1,j2], period)
            end
        end
        s2 /= n1*(n2 - 1)
    end

    # Compute the wrapped gradient (with mean slope subtracted) and the
    # curvature along each dimension.
    d1 = wrk1
    @inbounds for j2 in 1:n2
        @simd for j1 in 1:n1-1
            d1[j1,j2] = wrap(psi[j1+1,j2] - psi[j1,j2] - s1, period)
        end
        let j1 = n1
            d1[j1,j2] = 0
        end
        @simd for j1 in n1:-1:2
            d1[j1,j2] -= d1[j1-1,j2]
        end
    end
    d2 = wrk2
    @inbounds for j2 in 1:n2-1
        @simd for j1 in 1:n1
            d2[j1,j2] = wrap(psi[j1,j2+1] - psi[j1,j2] - s2, period)
        end
    end
    let j2 = n2
        @inbounds @simd for j1 in 1:n1
            d2[j1,j2] = 0
        end
    end
    @inbounds for j2 in n2:-1:2
        @simd for j1 in 1:n1
            d2[j1,j2] -= d2[j1,j2-1]
        end
    end

    # Compute the discrete Laplacian.
    rho = wrk1 # we re-use one of the workspaces to spare allocations
    @inbounds @simd for j1 in eachindex(wrk1, wrk2)
        rho[j1] = d1[j1] + d2[j1]
    end

    # Compute Q, the FFT filter implementing the discrete Laplacian.
    q = wrk2 # we re-use the other workspace to spare allocations
    q1 = Array{T}(undef, n1) # array of cosine along first dimension
    a1 = twopi(T)/n1
    @inbounds @simd for j1 in 1:n1
        q1[j1] = two_cos_minus_two(a1*(j1 - 1))
    end
    if n1 == n2
        @inbounds for j2 in 1:n2
            q2 = q1[j2]
            @simd for j1 in 1:n1
                q[j1,j2] = q1[j1] + q2
            end
        end
    else
        a2 = twopi(T)/n2
        @inbounds for j2 in 1:n2
            q2 = two_cos_minus_two(a2*(j2 - 1))
            @simd for j1 in 1:n1
                q[j1,j2] = q1[j1] + q2
            end
        end
    end
    q[1,1] = 1 # to avoid division by zero

    # Compute the solution.
    z = fft(rho)
    z[1,1] = 0 # should be exactly zero without rounding errors
    @inbounds @simd for j in eachindex(z, q)
        z[j] /= q[j]
    end
    ifft!(z)

    # Extract real part and restore the (centered) mean slopes.
    phi = wrk1
    c1 = T(n1 + 1)/2 # offset to have the slope centered
    c2 = T(n2 + 1)/2 # offset to have the slope centered
    @inbounds for j2 in 1:n2
        r2 = (j2 - c2)*s2
        @simd for j1 in 1:n1
            r1 = (j1 - c1)*s1
            phi[j1,j2] = z[j1,j2].re + (r1 + r2)
        end
    end

    # Fix the piston (modulo the period) in the computed solution.
    debiased && fix_bias!(phi, psi, period)

    return phi
end

"""
    fix_bias!(x, y, p=2π) -> x

fixes a constant bias modulo the period `p` in the array `x` which is an
unwrapped but biased version of `y`.  The operation is done in-place and `x` is
returned.

"""
function fix_bias!(x::AbstractArray{T,N},
                   y::AbstractArray{T,N},
                   period::T = twopi(T)) where {T<:AbstractFloat,N}
    a = twopi(T)/period
    sum_cos = zero(T)
    sum_sin = zero(T)
    @inbounds @simd for j in eachindex(x, y)
        sin_j, cos_j = sincos(a*(y[j]- x[j]))
        sum_sin += sin_j
        sum_cos += cos_j
    end
    bias = atan(sum_sin, sum_cos)/a
    @inbounds @simd for j in eachindex(x)
        x[j] += bias
    end
    return x
end

"""
    two_cos_minus_two(x) -> 2*cos(x) - 2

computes accurately `2*cos(x) - 2` given the floating-point value `x`.

"""
two_cos_minus_two(x::AbstractFloat) = -4*sin(x/2)^2

end # module
