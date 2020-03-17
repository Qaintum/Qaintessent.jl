using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(cg::CircuitGate)
    Qaintessent.matrix(cg) * Qaintessent.matrix(Base.adjoint(cg)) ≈ I
end

function isunitary(cb::CircuitBlock)
    Qaintessent.matrix(cb) * Qaintessent.matrix(Base.adjoint(cb)) ≈ I
end


@testset ExtendedTestSet "circuit gates" begin

    # Y acting on second wire
    cg = CircuitGate{1,3}((2,), Y)
    @test Qaintessent.matrix(cg) ≈ kron(kron(Matrix(I, 2, 2), Qaintessent.matrix(Y)), Matrix(I, 2, 2))
    @test isunitary(cg)

    # flip control and target
    cg = CircuitGate{2,2}((2, 1), controlled_not())
    @test Qaintessent.matrix(cg) ≈ [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]
    @test isunitary(cg)

    # first qubit as control and third qubit as target
    cg = controlled_circuit_gate(1, 3, HadamardGate(), 3)
    @test Qaintessent.matrix(cg) ≈ [
        Matrix(I, 4, 4) fill(0, 4, 2) fill(0, 4, 2);
        fill(0, 2, 4) Qaintessent.matrix(HadamardGate()) fill(0, 2, 2);
        fill(0, 2, 6) Qaintessent.matrix(HadamardGate())]
    @test isunitary(cg)
end


@testset ExtendedTestSet "simple circuit block" begin
    N = 1
    cb = CircuitBlock([
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        single_qubit_circuit_gate(1, Z, N),
        single_qubit_circuit_gate(1, Y, N),
    ])

    @test isunitary(cb)

    ψ = randn(ComplexF64, 2^N)
    @test Qaintessent.apply(cb, ψ) ≈ Qaintessent.matrix(cb)*ψ
    ψ = [1, 0, 0]
    @test_throws DimensionMismatch Qaintessent.apply(cb, ψ)

    # observable block with Pauli Z matrix
    ob = CircuitBlock([
        single_qubit_circuit_gate(1, Z, N),
        ])
    ψ = randn(ComplexF64, 2^N)
    ϕ1 = Qaintessent.measure(cb, ob, ψ)
    ϕ2 = Qaintessent.apply(cb, ψ)
    @test ϕ1 ≈ adjoint(ϕ2)* Qaintessent.apply(ob, ϕ2)

    # Arbitrary observerables. Non-physical, used to test math implementation.
    o = rand(2,2)
    ψ = randn(ComplexF64, 2^N)
    ϕ1 = Qaintessent.measure(cb, o, ψ)
    ϕ2 = Qaintessent.apply(cb, ψ)
    @test ϕ1 ≈ adjoint(ϕ2) * o * ϕ2
end


@testset ExtendedTestSet "three qubit QFT" begin
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

    @test cb[1] == single_qubit_circuit_gate(1, HadamardGate(), N)
    for (index, gate) in enumerate(cb)
        @test gate == cb[index]
    end

    @test Qaintessent.matrix(cb) ≈ [exp(2*π*1im*j*k/2^N)/sqrt(2^N) for j in 0:(2^N-1), k in 0:(2^N-1)]

    @test isunitary(cb)

    ψ = randn(ComplexF64, 2^N)
    @test Qaintessent.apply(cb, ψ) ≈ Qaintessent.matrix(cb)*ψ
end
