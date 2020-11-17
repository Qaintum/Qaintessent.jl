using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


@testset ExtendedTestSet "commute" begin

    # individual gates
    gates = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(2π*rand()), RyGate(2π*rand()), RzGate(2π*rand()), RotationGate(randn(3)), PhaseShiftGate(2π*rand()),
             SwapGate(), EntanglementXXGate(2π*rand()), EntanglementYYGate(2π*rand()), EntanglementZZGate(2π*rand()), controlled_not()]
    for g1 in gates
        for g2 in gates
            @test iscommuting(g1, g2) == iscommuting(Qaintessent.matrix(g1), Qaintessent.matrix(g2))
        end
    end

    # circuit gates
    N = 5
    ccg1 = begin
        iwperm = Tuple(randperm(N))
        # number of target and control wires
        nt = 1
        nc = rand(1:3)
        controlled_circuit_gate(iwperm[1:nt], iwperm[nt+1:nt+nc], HadamardGate(), N)
    end
    ccg2 = begin
        iwperm = Tuple(randperm(N))
        # number of target and control wires
        nt = 2
        nc = rand(1:3)
        controlled_circuit_gate(iwperm[1:nt], iwperm[nt+1:nt+nc], MatrixGate(Array(qr(randn(ComplexF64, 4, 4)).Q)), N)
    end
    cgates = [
        single_qubit_circuit_gate(rand(1:N), X, N),
        single_qubit_circuit_gate(rand(1:N), RzGate(2π*rand()), N),
        CircuitGate((2, 4), SwapGate(), N),
        ccg1,
        ccg2,
    ]
    for cg1 in cgates
        for cg2 in cgates
            @test iscommuting(cg1, cg2) == iscommuting(Qaintessent.matrix(cg1), Qaintessent.matrix(cg2))
        end
    end
end
