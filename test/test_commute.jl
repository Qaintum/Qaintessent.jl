using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "commute" begin

    # individual gates
    gates = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(2π*rand()), RyGate(2π*rand()), RzGate(2π*rand()), RotationGate(randn(3)), PhaseShiftGate(2π*rand()),
             SwapGate(), EntanglementXXGate(2π*rand()), EntanglementYYGate(2π*rand()), EntanglementZZGate(2π*rand()), controlled_not()]
    for g1 in gates
        for g2 in gates
            @test iscommuting(g1, g2) == iscommuting(sparse_matrix(g1), sparse_matrix(g2))
        end
    end

    # circuit gates
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
        circuit_gate(iwperm[1:nt], MatrixGate(Array(qr(randn(ComplexF64, 4, 4)).Q)), iwperm[nt+1:nt+nc],)
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
            @test iscommuting(cg1, cg2) == iscommuting(Qaintessent.sparse_matrix(cg1, N), Qaintessent.sparse_matrix(cg2, N))
        end
    end
end
