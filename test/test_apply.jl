using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


"""Checks that tailored apply gives same result as general apply"""

@testset ExtendedTestSet "apply" begin

    N = 6
    θ = 0.7*π
    ϕ = 0.4*π
    n = randn(3); n /= norm(n)
    ψ = rand(ComplexF64, 2^N)

    # one qubit gates
    for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
        i = rand(1:N)
        cg = single_qubit_circuit_gate(i, g, N)
        cga = CircuitGate{1,N,AbstractGate{1}}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}

        @test apply(cg, ψ) ≈ apply(cga, ψ)
    end

    # control gate
    for g in [X, Y, Z, RotationGate(θ,n), PhaseShiftGate(ϕ)]
        i = rand(1:N)
        j = rand([1:i-1; i+1:N])
        cg = controlled_circuit_gate(i, j, g, N)
        cga = CircuitGate{2,N,AbstractGate{2}}(cg.iwire, cg.gate)

        @test apply(cg, ψ) ≈ apply(cga, ψ)
    end

    # swap gate
    i = rand(1:N)
    j = rand([1:i-1; i+1:N])
    cg = CircuitGate((i,j), SwapGate(), N)
    cga = CircuitGate{2,N,AbstractGate{2}}(cg.iwire, cg.gate)

    @test apply(cg, ψ) ≈ apply(cga, ψ)
end


@testset ExtendedTestSet "density matrix apply" begin

    N = 5
    ψ = randn(ComplexF64, 2^N)
    ρ = density_from_statevector(ψ)

    for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(-1.1), RyGate(0.7), RzGate(0.4), RotationGate([-0.3, 0.1, 0.23]), PhaseShiftGate(0.9)]
        cg = CircuitGate((rand(1:N),), g, N)
        ψs = apply(cg, ψ)
        ρsref = density_from_statevector(ψs)
        ρs = apply(cg, ρ)
        @test ρs.v ≈ ρsref.v
    end
end
