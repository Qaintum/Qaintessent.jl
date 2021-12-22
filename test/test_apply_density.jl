using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent
using CUDA

##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "apply gates to density matrix" begin
    N = 5
    ψ = randn(ComplexQ, 2^N)
    ψ /= norm(ψ)
    # @code_warntype(density_from_statevector(ψ))
    

    @testset "density matrix apply single qubit gates" begin

        for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(-1.1), RyGate(0.7), RzGate(0.4), RotationGate([-0.3, 0.1, 0.23]), PhaseShiftGate(0.9), Qaintessent.Proj1Gate()]
            cg = circuit_gate(rand(1:N), g)
            ψs = apply(ψ, cg)
            
            ρ = density_from_statevector(ψ)
            ρsref = density_from_statevector(ψs)
            ρs1 = apply(ρ, cg)
            ρs2 = apply!(ρ, cg)
            # @code_warntype(apply(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v            

            ρ = density_from_statevector(ψ)
            ρsref = density_from_matrix(reshape(0.5 * (kron(conj(ψ), ψs) + kron(conj(ψs), ψ)), 2^N, 2^N))
            ρs1 = Qaintessent.apply_mixed_add(ρ, cg)            
            ρs2 = Qaintessent.apply_mixed_add!(ρ, cg)
            # @code_warntype(Qaintessent.apply_mixed_add(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v

            ρ = density_from_statevector(ψ)
            ρsref = density_from_matrix(reshape(0.5im * (kron(conj(ψ), ψs) - kron(conj(ψs), ψ)), 2^N, 2^N))
            ρs1 = Qaintessent.apply_mixed_sub(ρ, cg)
            ρs2 = Qaintessent.apply_mixed_sub!(ρ, cg)
            # @code_warntype(Qaintessent.apply_mixed_sub(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v
        end
    end

    @testset "density matrix apply two qubit gates" begin
        for g in [SwapGate(), EntanglementXXGate(-0.7), EntanglementYYGate(0.2*π), EntanglementZZGate(0.4)]
            i, j = randperm(N)[1:2]
            cg = CircuitGate((i, j), g)
            ψs = apply(ψ, cg)

            ρ = density_from_statevector(ψ)
            ρsref = density_from_statevector(ψs)
            ρs1 = apply(ρ, cg)
            ρs2 = apply!(ρ, cg)

            # @code_warntype(apply(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v     

            ρ = density_from_statevector(ψ)
            ρsref = density_from_matrix(reshape(0.5 * (kron(conj(ψ), ψs) + kron(conj(ψs), ψ)), 2^N, 2^N))
            ρs1 = Qaintessent.apply_mixed_add(ρ, cg)            
            ρs2 = Qaintessent.apply_mixed_add!(ρ, cg)
            # @code_warntype(Qaintessent.apply_mixed_add(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v

            ρ = density_from_statevector(ψ)
            ρsref = density_from_matrix(reshape(0.5im * (kron(conj(ψ), ψs) - kron(conj(ψs), ψ)), 2^N, 2^N))
            ρs1 = Qaintessent.apply_mixed_sub(ρ, cg)
            ρs2 = Qaintessent.apply_mixed_sub!(ρ, cg)
            # @code_warntype(Qaintessent.apply_mixed_sub(cg, ρ))

            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v
        end
    end

    @testset "density matrix apply controlled single qubit gate" begin
        # one target qubit
        for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(2π*rand()), RyGate(2π*rand()), RzGate(2π*rand()), RotationGate([2π*rand(), 2π*rand(), 2π*rand()]), PhaseShiftGate(2π*rand()), Qaintessent.Proj1Gate()]
            iwperm = randperm(N)
            cg = circuit_gate(iwperm[1], g, (iwperm[2], iwperm[3], iwperm[4]))

            ψs = apply(ψ, cg)
            ρsref = density_from_statevector(ψs)
            ρ = density_from_statevector(ψ)
            ρs1 = apply(ρ, cg)
            ρs2 = apply!(ρ, cg)
            # @code_warntype(apply(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v
        end
    end

    @testset "density matrix apply controlled two qubit gate" begin

        for g in [SwapGate(), EntanglementXXGate(2π*rand()), EntanglementYYGate(2π*rand()), EntanglementZZGate(2π*rand())]
            # two target qubits
            iwperm = randperm(N)
            cg = circuit_gate((iwperm[1], iwperm[2]), g, (iwperm[3], iwperm[4]))

            ψs = apply(ψ, cg)

            ρsref = density_from_statevector(ψs)
            ρ = density_from_statevector(ψ)
            ρs1 = apply(ρ, cg)
            ρs2 = apply!(ρ, cg)
            # @code_warntype(apply(cg, ρ))
            @test ρs1.v ≈ ρsref.v
            @test ρs2.v ≈ ρsref.v
        end
    end

    @testset "density matrix apply general unitary gate" begin
        # matrix gate (general unitary gate)
        cg = CircuitGate(NTuple{3,Int64}(randperm(N)[1:3]), MatrixGate(Matrix(qr(randn(ComplexQ, 8, 8)).Q)))
        ψs = apply(ψ, cg)
        ρ = density_from_statevector(ψ)
        ρsref = density_from_statevector(ψs)
        ρs1 = apply(ρ, cg)
        ρs2 = apply!(ρ, cg)
        # @code_warntype(apply(cg, ρ))
        @test ρs1.v ≈ ρsref.v
        @test ρs2.v ≈ ρsref.v
    end

    
    @testset "density matrix apply multiple unitary gates" begin
        iwperm = randperm(N)

        cgs = [circuit_gate((iwperm[1], iwperm[2]), EntanglementYYGate(2π*rand()), (iwperm[3], iwperm[4])),
            circuit_gate((iwperm[1], iwperm[2]), EntanglementZZGate(2π*rand()), (iwperm[3], iwperm[4])),
            circuit_gate((iwperm[1], iwperm[2]), EntanglementXXGate(2π*rand()), (iwperm[3], iwperm[4])),
        ]
        ψs = apply(ψ, cgs)

        ρsref = density_from_statevector(ψs)
        ρ = density_from_statevector(ψ)
        ρs1 = apply(ρ, cgs)
        ρs2 = apply!(ρ, cgs)
        # @code_warntype(apply(cg, ρ))
        @test ρs1.v ≈ ρsref.v
        @test ρs2.v ≈ ρsref.v
    end
end
