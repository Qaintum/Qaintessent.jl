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

    @testset "density matrix apply two qubit gates" begin
        for g in [SwapGate(), EntanglementXXGate(-0.7), EntanglementYYGate(0.2*π), EntanglementZZGate(0.4)]
            i, j = randperm(N)[1:2]
            cg = CircuitGate((i, j), g)
            ψs = apply(cg, ψ)

            ρsref = density_from_statevector(ψs)
            ρs = apply(cg, ρ)
            #@code_warntype(apply(cg, ρ))
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
    end

    @testset "density matrix apply controlled gate" begin
        # one target qubit
        iwperm = randperm(N)
        cg = circuit_gate(iwperm[1], RotationGate(rand(3) .- 0.5), (iwperm[2], iwperm[3], iwperm[4]))

        ψs = apply(cg, ψ)
        ρsref = density_from_statevector(ψs)

        ρs = apply(cg, ρ)
        # @code_warntype(apply(cg, ρ))
        @test ρs.v ≈ ρsref.v

        # two target qubits
        iwperm = randperm(N)
        cg = circuit_gate((iwperm[1], iwperm[2]), EntanglementYYGate(2π*rand()), (iwperm[3], iwperm[4]))

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
        # @code_warntype(apply(cg, ρ))
        @test ρs.v ≈ ρsref.v
    end
end
