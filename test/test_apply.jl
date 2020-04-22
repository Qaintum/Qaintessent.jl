using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent

""" Checks that tailored apply gives same result as general apply """

@testset ExtendedTestSet "apply" begin

    N = 6
    θ = 0.7*π
    ϕ = 0.4*π
    n = randn(3); n /= norm(n)
    ψ = rand(ComplexF64,2^N)

    # one qubit gates
    for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
        i = rand(1:N)
        U = single_qubit_circuit_gate(i, g, N)
        U_AbsGate = CircuitGate{1,N,AbstractGate}(U.iwire,U.gate) # generate same gate with type AbstractGate

        @test apply(U, ψ) ≈ apply(U_AbsGate, ψ)
    end

    # control gate
    for g in [X, Y, Z, RotationGate(θ,n), PhaseShiftGate(ϕ)]
        i = rand(1:N)
        j = rand([1:i-1; i+1:N])
        U = controlled_circuit_gate(i, j, g, N)
        U_AbsGate = CircuitGate{2,N,AbstractGate}(U.iwire,U.gate)

        @test apply(U, ψ) ≈ apply(U_AbsGate, ψ)
    end

    # swap gate
    i = rand(1:N)
    j = rand([1:i-1; i+1:N])
    U = CircuitGate{2,N}((i,j),SwapGate())
    U_AbsGate = CircuitGate{2,N,AbstractGate}(U.iwire,U.gate)

    @test apply(U, ψ) ≈ apply(U_AbsGate, ψ)


end
