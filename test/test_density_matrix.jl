using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "density matrix" begin

    pauli = [
        Matrix{Float64}(I, 2, 2),
        Qaintessent.matrix(X),
        Qaintessent.matrix(Y),
        Qaintessent.matrix(Z),
    ]

    @test Qaintessent.matrix(pauli_group_matrix("XIZY")) ≈
          kron(pauli[3], pauli[4], pauli[1], pauli[2])

    N = 4
    ψ = randn(ComplexF64, 2^N)
    # matrix representation of density matrix |ψ⟩⟨ψ|
    ρmat = reshape(kron(conj(ψ), ψ), 2^N, 2^N)
    # density matrix represented with respect to Pauli basis
    ρ = density_from_statevector(ψ)
    @test Qaintessent.matrix(ρ) ≈ ρmat
    # compute Pauli representation from matrix form
    ρalt = density_from_matrix(ρmat)
    @test ρalt.v ≈ ρ.v
end

@testset ExtendedTestSet "density matrix helper functions" begin
    N = 4
    ψ1 = randn(ComplexF64, 2^N)
    ψ2 = randn(ComplexF64, 2^N)
    ρ1 = density_from_statevector(ψ1)
    ρ2 = density_from_statevector(ψ2)

    @test dot(ρ1, ρ2) ≈ dot(ρ1.v, ρ2.v)
end