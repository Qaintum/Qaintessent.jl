using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "density matrix" begin

    pauli = [
        Matrix{FloatQ}(I, 2, 2),
        Qaintessent.matrix(X),
        Qaintessent.matrix(Y),
        Qaintessent.matrix(Z),
    ]
    @testset "density matrix pauli group matrix" begin
        @test Qaintessent.matrix(pauli_group_matrix("XIZY")) ≈
            kron(pauli[3], pauli[4], pauli[1], pauli[2])
    end

    @testset "density matrix constructor" begin
        N = 3

        v = randn(FloatQ, 4^N)
        d1 = DensityMatrix(v)
        d2 = DensityMatrix(v, N)
        @test d1.v ≈ d2.v
        @test d1.N == d2.N
    end

    @testset "density matrix constructor helper function" begin
        N = 4
        ψ = randn(ComplexQ, 2^N)
        # matrix representation of density matrix |ψ⟩⟨ψ|
        ρmat = reshape(kron(conj(ψ), ψ), 2^N, 2^N)
        # density matrix represented with respect to Pauli basis
        ρ = density_from_statevector(ψ)
        @test Qaintessent.matrix(ρ) ≈ ρmat
        # compute Pauli representation from matrix form
        ρalt = density_from_matrix(ρmat)
        @test ρalt.v ≈ ρ.v
    end

end
