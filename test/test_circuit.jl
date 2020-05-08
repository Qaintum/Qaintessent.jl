using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(cg::CircuitGate)
    Qaintessent.matrix(cg) * Qaintessent.matrix(Base.adjoint(cg)) ≈ I
end

function isunitary(cgc::CircuitGateChain)
    Qaintessent.matrix(cgc) * Qaintessent.matrix(Base.adjoint(cgc)) ≈ I
end

@testset ExtendedTestSet "moments" begin
    N = 3
    @test_throws ErrorException cgc_ref = CircuitGateChain{N}([
        Moment{N}(
        [single_qubit_circuit_gate(2, X, N),
        single_qubit_circuit_gate(2, HadamardGate(), N),
        single_qubit_circuit_gate(3, HadamardGate(), N)]
        ),
        Moment{N}(
        [single_qubit_circuit_gate(1, Z, N),
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(1, Y, N)]
        ),
    ])
end

@testset ExtendedTestSet "circuit gates" begin
    θ = 0.7*π
    ϕ = 0.4*π
    n = randn(3); n /= norm(n)

    # single qubit gates
    for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
        cg = CircuitGate((2,), g, 3)
        @test LinearAlgebra.ishermitian(cg) == (Qaintessent.matrix(cg) == Qaintessent.matrix(adjoint(cg)))
    end

    # two qubit gates
    for g in [controlled_not()]
        cg = CircuitGate((2, 3), g, 3)
        @test LinearAlgebra.ishermitian(cg) == (Qaintessent.matrix(cg) == Qaintessent.matrix(adjoint(cg)))
    end

    # Y acting on second wire
    cg = CircuitGate((2,), Y, 3)
    @test Qaintessent.matrix(cg) ≈ kron(kron(Matrix(I, 2, 2), Qaintessent.matrix(Y)), Matrix(I, 2, 2))
    @test isunitary(cg)

    # flip control and target
    cg = CircuitGate((2, 1), controlled_not(), 2)
    @test Qaintessent.matrix(cg) ≈ [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]
    @test isunitary(cg)

    # first qubit as control and third qubit as target
    cg = controlled_circuit_gate(1, 3, HadamardGate(), 3)
    @test Qaintessent.matrix(cg) ≈ [
        Matrix(I, 4, 4) fill(0, 4, 2) fill(0, 4, 2);
        fill(0, 2, 4) Qaintessent.matrix(HadamardGate()) fill(0, 2, 2);
        fill(0, 2, 6) Qaintessent.matrix(HadamardGate())]
    @test isunitary(cg)

    cg_test = controlled_circuit_gate(1, 3, HadamardGate(), 3)
    @test cg_test == cg

    cg = single_qubit_circuit_gate(1, RxGate(0.2), 1)
    cg_test = single_qubit_circuit_gate(1, RxGate(0.2), 1)

    @test cg_test ≈ cg

    cg.gate.θ[1] = 0.5
    @test !(cg ≈ cg_test)
end

@testset ExtendedTestSet "circuit gate chain" begin
    N = 1
    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        single_qubit_circuit_gate(1, Z, N),
        single_qubit_circuit_gate(1, Y, N),
    ])

    cgc1 = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
    ])

    cgc2 = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, Z, N),
        single_qubit_circuit_gate(1, Y, N),
    ])

    cgc_test = cgc1 * cgc2

    for (gate1, gate2) in zip(cgc_test, cgc_ref)
        @test gate1 ≈ gate2
    end
end

@testset ExtendedTestSet "circuit" begin
    N = 1
    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        single_qubit_circuit_gate(1, Z, N),
        single_qubit_circuit_gate(1, Y, N),
    ])
    meas = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N)])
    c = Circuit(cgc, meas)
    for i in 1:length(cgc)
        @test c[i] == cgc[i]
    end
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.cgc, ψ)
end

@testset ExtendedTestSet "reduced density matrix" begin
    N = 4
    ψ = randn(ComplexF64, 2^N)
    χ = randn(ComplexF64, 2^N)
    # full density matrix |ψ⟩⟨χ|
    ρ = reshape(kron(conj(χ), ψ), 2^N, 2^N)

    id = Matrix{ComplexF64}(I, 2, 2)
    A = randn(ComplexF64, 2, 2)
    B = randn(ComplexF64, 2, 2)

    @test sum(kron(A, B) .* rdm(N, (4, 2), ψ, χ)) ≈ sum(kron(id, B, id, A) .* ρ)
end


@testset ExtendedTestSet "simple circuit chain" begin
    N = 1
    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        single_qubit_circuit_gate(1, Z, N),
        single_qubit_circuit_gate(1, Y, N),
    ])

    @test isunitary(cgc)

    ψ = randn(ComplexF64, 2^N)
    @test apply(cgc, ψ) ≈ Qaintessent.matrix(cgc)*ψ
end


@testset ExtendedTestSet "three qubit QFT" begin
    # three qubit quantum Fourier transform
    N = 3
    cgc = CircuitGateChain{N}([
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
    @test cgc[1] ≈ Moment{N}(single_qubit_circuit_gate(1, HadamardGate(), N))
    for (index, moment) in enumerate(cgc)
        @test moment == cgc[index]
    end

    @test Qaintessent.matrix(cgc) ≈ [exp(2*π*1im*j*k/2^N)/sqrt(2^N) for j in 0:(2^N-1), k in 0:(2^N-1)]

    @test isunitary(cgc)

    ψ = randn(ComplexF64, 2^N)
    @test apply(cgc, ψ) ≈ Qaintessent.matrix(cgc)*ψ
end


@testset ExtendedTestSet "quantum circuit" begin
    N = 3

    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(2, Y, N),
        controlled_circuit_gate((1, 3), 2, SGate(), N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(3, RxGate(1.5π), N),
    ])

    Xmat = Qaintessent.matrix(X)
    Ymat = Qaintessent.matrix(Y)
    Zmat = Qaintessent.matrix(Z)
    mop1 = kron(Matrix{Float64}(I, 2, 2), kron(Xmat, Xmat))
    mop2 = kron(Matrix{Float64}(I, 2, 2), kron(Zmat, Ymat))
    mop3 = kron(Zmat, Matrix{Float64}(I, 4, 4))
    meas = MeasurementOps{N}([mop1, mop2, mop3])

    c = Circuit(cgc, meas)

    ψ = randn(ComplexF64, 2^N)

    ψs = apply(cgc, ψ)
    @test [dot(ψs, m*ψs) for m in [mop1, mop2, mop3]] ≈ apply(c, ψ)

end
