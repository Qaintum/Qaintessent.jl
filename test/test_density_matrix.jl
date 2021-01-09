using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


@testset ExtendedTestSet "density matrix" begin

    pauli = [
        Matrix{Float64}(I, 2, 2),
        Qaintessent.sparse_matrix(X),
        Qaintessent.sparse_matrix(Y),
        Qaintessent.sparse_matrix(Z),
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
end
