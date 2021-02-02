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

    @testset "density matrix apply single qubit gates" begin

        for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(-1.1), RyGate(0.7), RzGate(0.4), RotationGate([-0.3, 0.1, 0.23]), PhaseShiftGate(0.9), Qaintessent.Proj1Gate()]
            cg = circuit_gate(rand(1:N), g)
            ψs = apply(cg, ψ)

            ρsref = density_from_statevector(ψs)
            ρs = apply(cg, ρ)
            # @code_warntype(apply(cg, ρ))
            @test ρs.v ≈ ρsref.v

            # generate same gate with type AbstractGate
            cga = CircuitGate{1,AbstractGate}(cg.iwire, cg.gate)
            ρsa = apply(cga, ρ)
            # @code_warntype(apply(cga, ρ))
            @test ρs.v ≈ ρsa.v

            ρsref = density_from_matrix(reshape(0.5 * (kron(conj(ψ), ψs) + kron(conj(ψs), ψ)), 2^N, 2^N))
            ρs = Qaintessent.apply_mixed_add(cg, ρ)
            # @code_warntype(Qaintessent.apply_mixed_add(cg, ρ))
            @test ρs.v ≈ ρsref.v

            ρsref = density_from_matrix(reshape(0.5im * (kron(conj(ψ), ψs) - kron(conj(ψs), ψ)), 2^N, 2^N))
            ρs = Qaintessent.apply_mixed_sub(cg, ρ)
            # @code_warntype(Qaintessent.apply_mixed_sub(cg, ρ))
            @test ρs.v ≈ ρsref.v
        end
    end

    @testset "density matrix apply swap gate" begin
        # swap gate
        i, j = randperm(N)[1:2]
        cg = CircuitGate((i, j), SwapGate())
        ψs = apply(cg, ψ)

        ρsref = density_from_statevector(ψs)
        ρs = apply(cg, ρ)
        # @code_warntype(apply(cg, ρ))
        @test ρs.v ≈ ρsref.v

        ρsref = density_from_matrix(reshape(0.5 * (kron(conj(ψ), ψs) + kron(conj(ψs), ψ)), 2^N, 2^N))
        ρs = Qaintessent.apply_mixed_add(cg, ρ)
        # @code_warntype(Qaintessent.apply_mixed_add(cg, ρ))
        @test ρs.v ≈ ρsref.v

        ρsref = density_from_matrix(reshape(0.5im * (kron(conj(ψ), ψs) - kron(conj(ψs), ψ)), 2^N, 2^N))
        ρs = Qaintessent.apply_mixed_sub(cg, ρ)
        # @code_warntype(Qaintessent.apply_mixed_sub(cg, ρ))
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
        # @code_warntype(apply(cg, ρ))
        @test ρs.v ≈ ρsref.v
    end

    @testset "density matrix apply general unitary gate" begin
        # matrix gate (general unitary gate)
        cg = CircuitGate(NTuple{3,Int64}(randperm(N)[1:3]), MatrixGate(Matrix(qr(randn(ComplexF64, 8, 8)).Q)))
        ψs = apply(cg, ψ)
        ρsref = density_from_statevector(ψs)
        ρs = apply(cg, ρ)
        @test ρs.v ≈ ρsref.v
    end
end
