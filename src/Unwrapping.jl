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
    unwrap(psi, p=2π) -> phi

yields the unwrapped phase `phi` that best fits, in a least squares sense, the
wrapped phase `psi`.  Optional argument `p` is the period of the wrapping.
The result has an undetermined constant offset.

The method is based on algorithm *"Robust two-dimensional weighted and
unweighted phase unwrapping that uses fast transforms and iterative methods"*
by Dennis C. Ghiglia and Louis A. Romero (J. Opt. Soc. Am. A, Vol. 11,
pp. 107-117, 1994).

"""
function unwrap(psi::AbstractMatrix{T},
                period::Real; kwds...) where {T<:AbstractFloat}
    unwrap(psi, T(period); kwds...)
end

function unwrap(psi::AbstractMatrix{T},
                period::T = twopi(T)) where {T<:AbstractFloat}
    Base.has_offset_axes(psi) && error("input array has non-standard axes")
    n1, n2 = size(psi)

    # We need 2 workspace arrays.
    wrk1 = fill!(Array{T}(undef, n1, n2),0)
    wrk2 = fill!(Array{T}(undef, n1, n2),0)

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
    rho = wrk1 # we re-use one of the workspace to spare allocations
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
        q1[j1] = 2*(cos(a1*(j1 - 1)) - 1)
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
            q2 = 2*(cos(a2*(j2 - 1)) - 1)
            @simd for j1 in 1:n1
                q[j1,j2] = q1[j1] + q2
            end
        end
    end
    q[1,1] = 1 # to avoid division by zero

    # Compute the solution (Julia's dct is of type II).
    dct!(rho)
    phi = wrk1 # we re-use one of the workspaces to spare allocations
    @inbounds @simd for j in eachindex(wrk1, wrk2)
        phi[j] = rho[j]/q[j]
    end
    idct!(phi)
    return phi
end

end # module
