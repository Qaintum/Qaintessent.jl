using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


@testset ExtendedTestSet "apply gates to density matrix" begin
    N = 5
    ψ = randn(ComplexF64, 2^N)
    ψ /= norm(ψ)
    ρ = density_from_statevector(ψ)

    @testset "density matrix apply basic gates" begin
        # single qubit gates
        for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(-1.1), RyGate(0.7), RzGate(0.4), RotationGate([-0.3, 0.1, 0.23]), PhaseShiftGate(0.9)]
            cg = CircuitGate((rand(1:N),), g)
            ψs = apply(cg, ψ)
            ρsref = density_from_statevector(ψs)
            ρs = apply(cg, ρ)
            @test ρs.v ≈ ρsref.v

            # generate same gate with type AbstractGate{1}
            cga = CircuitGate{1,AbstractGate}(cg.iwire, cg.gate)
            ρsa = apply(cga, ρ)
            @test ρs.v ≈ ρsa.v
        end
    end

    @testset "density matrix apply swap gate" begin
        # swap gate
        i = rand(1:N)
        j = rand([1:i-1; i+1:N])
        cg = CircuitGate((i, j), SwapGate())
        ψs = apply(cg, ψ)
        ρsref = density_from_statevector(ψs)
        ρs = apply(cg, ρ)
        @test ρs.v ≈ ρsref.v
    end

    @testset "density matrix apply controlled gate" begin
        # controlled gate
        iwperm = Tuple(randperm(N))
        # number of control and target wires
        nt = rand(1:2)
        nc = rand(1:3)
        cg = circuit_gate(iwperm[1:nt], nt == 1 ? RotationGate(rand(3) .- 0.5) : MatrixGate(Array(qr(randn(ComplexF64, 4, 4)).Q)), iwperm[nt+1:nt+nc])
        ψs = apply(cg, ψ)
        ρsref = density_from_statevector(ψs)
        ρs = apply(cg, ρ)
        @test ρs.v ≈ ρsref.v
    end

    @testset "density matrix apply general unitary gate" begin
        # matrix gate (general unitary gate)
        iwperm = Tuple(randperm(N))
        cg = CircuitGate(iwperm[1:3], MatrixGate(Array(qr(randn(ComplexF64, 8, 8)).Q)))
        ψs = apply(cg, ψ)
        ρsref = density_from_statevector(ψs)
        ρs = apply(cg, ρ)
        @test ρs.v ≈ ρsref.v
    end
end
