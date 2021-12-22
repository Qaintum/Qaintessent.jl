using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "commute" begin

    gates = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(2π*rand()), RyGate(2π*rand()), RzGate(2π*rand()), RotationGate(randn(3)), PhaseShiftGate(2π*rand()),
    SwapGate(), EntanglementXXGate(2π*rand()), EntanglementYYGate(2π*rand()), EntanglementZZGate(2π*rand()), controlled_not()]
    # individual gates
    @testset "commute basic gates" begin  
        for g1 in gates
            for g2 in gates
                @test iscommuting(g1, g2) == iscommuting(sparse_matrix(g1), sparse_matrix(g2))
            end
        end
    end
    # circuit gates
    @testset "commute circuit gates" begin
        N = 5
        ccg1 = begin
            iwperm = Tuple(randperm(N))
            # number of target and control wires
            nt = 1
            nc = rand(1:3)
            circuit_gate(iwperm[1:nt], HadamardGate(), iwperm[nt+1:nt+nc])
        end
        ccg2 = begin
            iwperm = Tuple(randperm(N))
            # number of target and control wires
            nt = 2
            nc = rand(1:3)
            circuit_gate(iwperm[1:nt], MatrixGate(Array(qr(randn(ComplexQ, 4, 4)).Q)), iwperm[nt+1:nt+nc],)
        end
        cgates = [
            circuit_gate(rand(1:N), X),
            circuit_gate(rand(1:N), RzGate(2π*rand())),
            CircuitGate((2, 4), SwapGate()),
            ccg1,
            ccg2,
        ]
        for cg1 in cgates
            for cg2 in cgates
                @test iscommuting(cg1, cg2) == iscommuting(sparse_matrix(cg1, N), sparse_matrix(cg2, N))
            end
        end
    end

    @testset "commute controlled circuit gates" begin
        N = 6
        iwire1, iwire2, iwire3, iwire4 = sample(1:N, 4, replace=false)
        cgates = [
            circuit_gate(iwire1, X, iwire2),
            circuit_gate(iwire3, RzGate(2π*rand()), iwire4),
            circuit_gate((iwire1, iwire4), SwapGate(), iwire2),
            circuit_gate(iwire2, Z, iwire3)
        ]
        for cg1 in cgates
            for cg2 in cgates
                @test iscommuting(cg1, cg2) == iscommuting(sparse_matrix(cg1, N), sparse_matrix(cg2, N))
            end
        end

    end

    @testset "commute measurement operators" begin
        N = 4
        q1 = Array(qr(randn(ComplexQ, 4, 4)).Q)
        q2 = Array(qr(randn(ComplexQ, 4, 4)).Q)
        q3 = Array(qr(randn(ComplexQ, 4, 4)).Q)
        mops = [
            mop(X, rand(1:N)),
            mop(Y, rand(1:N)),
            mop(SwapGate(), (2, 4)),
            mop(q1+q1', (3, 1)),
            mop(q2+q2', (2, 4)),
            mop(q3+q3', (1, 2))
        ]

        for mop1 in mops
            for mop2 in mops
                @test iscommuting(mop1, mop2) == iscommuting(sparse_matrix(mop1, N), sparse_matrix(mop2, N))
            end
        end
    end
end
