module UnwrappingTests

using Unwrapping, Test, Statistics
using Unwrapping: unwrap, wrap, arc

function make_screen(T::Type{<:AbstractFloat}, dims::Dims{2})
    A = Array{T}(undef, dims)
    n1, n2 = dims
    p = Unwrapping.twopi(T) # period
    c1 = T(n1+1)/2
    c2 = T(n2+1)/2
    # keep |ai|*ni + |bi| < p to avoid excessive wrapping
    a1 = p/(3*n1)
    a2 = p/(4*n2)
    b1 = T(-0.4)
    b2 = T(+0.3)
    for i2 in 1:n2
        x2 = i2 - c2
        for i1 in 1:n1
            x1 = i1 - c1
            A[i1,i2] = (a1*x1 + b1)*x1 + (a2*x2 + b2)*x2
        end
    end
    return A
end

dims = (384, 288)

@testset "Original Ghiglia & Romero method (T = $T)" for T in (Float32, Float64)
    tol = eps(T)*8
    A = make_screen(T, dims)
    Aw = arc.(A)
    A1 = unwrap(Aw; debiased=false)
    A2 = unwrap(Aw; debiased=true)
    # std() method is insensitive to a bias
    @test std(A - A1) ≤ tol*std(A - Aw)
    @test std(A - A2) ≤ tol*std(A - Aw)
    # A2 is debiased (modulo 2π)
    @test abs(arc(mean(A - A2))) ≤ sqrt(tol)
end

end # module
