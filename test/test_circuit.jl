using Test
using LinearAlgebra

using Qaintessent


@testset "circuit gates" begin

    # Y acting on second wire
    cg = CircuitGate{1,3}((2,), Y)
    @test Qaintessent.matrix(cg) ≈ kron(kron(Matrix(I, 2, 2), Qaintessent.matrix(Y)), Matrix(I, 2, 2))

    # flip control and target
    cg = CircuitGate{2,2}((2, 1), controlled_not())
    @test Qaintessent.matrix(cg) ≈ [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]

    # first qubit as control and third qubit as target
    cg = controlled_circuit_gate(1, 3, HadamardGate(), 3)
    @test Qaintessent.matrix(cg) ≈ [
        Matrix(I, 4, 4) fill(0, 4, 2) fill(0, 4, 2);
        fill(0, 2, 4) Qaintessent.matrix(HadamardGate()) fill(0, 2, 2);
        fill(0, 2, 6) Qaintessent.matrix(HadamardGate())]

end


@testset "circuit block" begin

    # three qubit quantum Fourier transform
    N = 3
    cb = CircuitBlock([
        # first Hadamard gate
        single_qubit_circuit_gate(1, HadamardGate(), N),
        # first controlled-S gate
        controlled_circuit_gate(2, 1, SGate(), N),
        # controlled-T gate
        controlled_circuit_gate(3, 1, TGate(), N),
        # second Hadamard gate
        single_qubit_circuit_gate(2, HadamardGate(), N),
        # second controlled-S gate
        controlled_circuit_gate(3, 2, SGate(), N),
        # third Hadamard gate
        single_qubit_circuit_gate(3, HadamardGate(), N),
        # final swap gate
        two_qubit_circuit_gate(1, 3, SwapGate(), N),
    ])

    @test Qaintessent.matrix(cb) ≈ [exp(2*π*1im*j*k/8)/sqrt(8) for j in 0:7, k in 0:7]

end
