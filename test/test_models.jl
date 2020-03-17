using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(cb::CircuitBlock)
    Qaintessent.matrix(cb) * Qaintessent.matrix(Base.adjoint(cb)) ≈ I
end


@testset ExtendedTestSet "quantum Fourier transform" begin
    for N in 1:4
        c = qft_circuit(N)

        Fmat = [exp(2*π*im*j*k/2^N)/sqrt(2^N) for j in 0:(2^N-1), k in 0:(2^N-1)]
        @test Qaintessent.matrix(c) ≈ Fmat

        @test isunitary(c)
    end
end
