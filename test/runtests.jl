module UnwrappingTests

using Unwrapping, Test, Statistics, GZip
using Unwrapping: unwrap, fft_unwrap, wrap, arc

const TESTDIR = @__DIR__

function read_test_data(filename::AbstractString)
    name = basename(filename)
    gzipped = endswith(name, ".gz")
    last = length(name)
    if gzipped
        last -= 3
    end
    i3 = findprev(isequal('.'), name, last)
    (i3 === nothing || SubString(name, i3:i3+3) != ".dat") &&
        error("invalid file extension")
    i2 = findprev(isequal('x'), name, i3-1)
    i2 === nothing && error("invalid file name (missing 'x' in dimensions)")
    i1 = findprev(isequal('-'), name, i2-1)
    i1 === nothing && error("invalid file name (missing dimensions)")
    n1 = tryparse(Int, SubString(name, i1+1:i2-1); base=10)
    (n1 === nothing || n1 < 1) &&
        error("invalid file name (bad first dimension)")
    n2 = tryparse(Int, SubString(name, i2+1:i3-1); base=10)
    (n2 === nothing || n2 < 1) &&
        error("invalid file name (bad second dimension)")
    A = Array{UInt8}(undef, n1, n2)
    if gzipped
        GZip.open(filename, "r") do io
            read!(io, A)
        end
    else
        open(filename, "r") do io
            read!(io, A)
        end
    end
    return A
end

function make_screen(T::Type{<:AbstractFloat}, dims::Dims{2})
    A = Array{T}(undef, dims)
    n1, n2 = dims
    p = Unwrapping.twopi(T) # period
    c1 = T(n1+1)/2
    c2 = T(n2+1)/2
    # keep |ai|*ni + |bi| << p to avoid excessive wrapping
    a1 = p/(7*n1)
    a2 = p/(9*n2)
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

gauss = read_test_data(joinpath(TESTDIR, "gaussian-434x342.dat.gz"))
peaks = read_test_data(joinpath(TESTDIR, "peaks-434x342.dat.gz"))

#dims = (384, 288)

@testset "DCT and FFT-based methods (T = $T)" for T in (Float32, Float64)
    dct_tol = 200*eps(T)
    fft_tol = 0.005
    #A = make_screen(T, dims)
    for A in (T(0.2)*gauss, T(0.2)*peaks)
        Aw = arc.(A)
        A1 = unwrap(Aw; debiased=false)
        A2 = unwrap(Aw; debiased=true)
        B1 = fft_unwrap(Aw; debiased=false)
        B2 = fft_unwrap(Aw; debiased=true)

        # std() method is insensitive to a bias
        @test std(A - A1) ≤ dct_tol*std(A - Aw)
        @test std(A - A2) ≤ dct_tol*std(A - Aw)
        @test std(A - B1) ≤ fft_tol*std(A - Aw)
        @test std(A - B2) ≤ fft_tol*std(A - Aw)
        # A2 and B2 are debiased (modulo 2π)
        @test abs(arc(mean(A - A2))) ≤ 2e-3
        @test abs(arc(mean(A - B2))) ≤ 2e-3
    end
end

end # module
