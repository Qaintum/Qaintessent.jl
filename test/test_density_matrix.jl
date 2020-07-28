using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


@testset ExtendedTestSet "density matrix" begin

    pauli = [
        Matrix{Float64}(I, 2, 2),
        Qaintessent.matrix(X),
        Qaintessent.matrix(Y),
        Qaintessent.matrix(Z),
    ]

    @test Qaintessent.matrix(pauli_group_matrix("XIZY")) ≈
          kron(pauli[2], pauli[1], pauli[4], pauli[3])

    N = 4
    ψ = randn(ComplexF64, 2^N)
    # matrix representation of density matrix |ψ⟩⟨ψ|
    ρmat = reshape(kron(conj(ψ), ψ), 2^N, 2^N)
    # density matrix represented with respect to Pauli basis
    ρ = density_from_statevector(ψ)
    @test Qaintessent.matrix(ρ) ≈ ρmat
end
