using Test
using LinearAlgebra

include("../src/Qaintessent.jl")
using .Qaintessent


@testset "quantum circuits" begin

    # Y acting on second wire
    cg = CircuitGate{1,3}((2,), Y)
    @test Qaintessent.matrix(cg) ≈ kron(kron(Matrix(I, 2, 2), Qaintessent.matrix(Y)), Matrix(I, 2, 2))

    # flip control and target
    cg = CircuitGate{2,2}((2, 1), controlled_not())
    @test Qaintessent.matrix(cg) ≈ [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]

    cg = CircuitGate{2,3}((1, 3), ControlledGate{1,1,2}(HadamardGate()))
    @test Qaintessent.matrix(cg) ≈ [
        Matrix(I, 4, 4) fill(0, 4, 2) fill(0, 4, 2);
        fill(0, 2, 4) Qaintessent.matrix(HadamardGate()) fill(0, 2, 2);
        fill(0, 2, 6) Qaintessent.matrix(HadamardGate())]

end
