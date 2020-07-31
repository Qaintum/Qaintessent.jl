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
    @test cg_test ≈ cg

    cg = single_qubit_circuit_gate(1, RxGate(0.2), 1)
    cg_test = single_qubit_circuit_gate(1, RxGate(0.2), 1)

    @test cg_test ≈ cg

    cg.gate.θ[1] = 0.5
    @test !(cg ≈ cg_test)
end

@testset ExtendedTestSet "circuit gate exceptions" begin
    H = HadamardGate()
    S = SwapGate()

    N = 1
    @test_throws ErrorException("Number of gate wires cannot be larger than total number of wires.") CircuitGate{2,N,SwapGate}((1, 2), S)

    N = 2
    @test_throws ErrorException("Need at least one wire to act on.") CircuitGate{0,N,SwapGate}((), S)
    @test_throws ErrorException("Wire indices must be unique.") CircuitGate{2,N,SwapGate}((1, 1), S)
    @test_throws ErrorException("Wire index cannot be smaller than 1.") CircuitGate{2,N,SwapGate}((1, -1), S)
    @test_throws ErrorException("Wire index cannot be larger than total number of wires.") CircuitGate{2,N,SwapGate}((1, 4), S)
    @test_throws ErrorException("Gate type must be a subtype of AbstractGate{1}.") CircuitGate{1,N,SwapGate}((1,), S)

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

    @test all(cgc_test .≈ cgc_ref)

    ψ = randn(ComplexF64, 2^N)
    ψ_ref = deepcopy(ψ)
    for c in cgc_ref
        ψ_ref = apply(c, ψ_ref)
    end

    @test ψ_ref ≈ apply(cgc_test, ψ)
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

    ψs = apply(cgc, ψ)
    @test [dot(ψs, m*ψs) for m in meas.mops] ≈ apply(c, ψ)

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


@testset ExtendedTestSet "CRegister" begin
    N = 5

    c1 = creg(2)
    c2 = creg(4)

    cgc = CircuitGateChain{N}([c1, c2])

    Qaintessent.set_creg!(c1, 2)
    Qaintessent.set_creg!(c2, 1, true)

    gates = [
        single_qubit_circuit_gate(2, Y, N),
        controlled_circuit_gate(reg_check(c2, 3), (2,), SGate(), N),
        controlled_circuit_gate((c2[2], 3), 2, YGate(), N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(5, RxGate(1.5π), N),
        single_qubit_circuit_gate(3, RyGate(1.5π), N),
        controlled_circuit_gate(c2[1], 2, ZGate(), N),
    ]
    cgc(gates)

    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(2, Y, N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(5, RxGate(1.5π), N),
        single_qubit_circuit_gate(3, RyGate(1.5π), N),
        single_qubit_circuit_gate(2, ZGate(), N)
    ])
    ψ = randn(ComplexF64, 2^N)
    @test apply(cgc_ref, ψ) ≈ apply(cgc, ψ)

    c3 = qreg(4)

end

@testset ExtendedTestSet "test add_creg!" begin
    N = 7
    c1 = creg(2)
    c2 = creg(4)
    c3 = creg(3)

    cgc = CircuitGateChain{N}()

    add_creg!(cgc, c1)
    add_creg!(cgc, c2)

    @test_throws ErrorException("CRegister is already used in another CircuitGateChain object") add_creg!(cgc, c2)
    @test_throws ErrorException("Register object has yet to be used in a CircuitGateChain object") gates = [controlled_circuit_gate(c3[1], 2, SGate(), N),]

    gates = [
        single_qubit_circuit_gate(2, Y, N),
        controlled_circuit_gate(reg_check(c2,7), 2, SGate(), N),
        controlled_circuit_gate((c2[4], 3), 2, YGate(), N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(5, RxGate(1.5π), N),
        single_qubit_circuit_gate(3, RyGate(1.5π), N),
        controlled_circuit_gate((c2[1]), 2, ZGate(), N),
    ]
    cgc(gates)

    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(2, SGate(), N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(5, RxGate(1.5π), N),
        single_qubit_circuit_gate(3, RyGate(1.5π), N),
        single_qubit_circuit_gate(2, ZGate(), N)
    ])

    ψ = randn(ComplexF64, 2^N)
    @test !(apply(cgc_ref, ψ) ≈ apply(cgc, ψ))
    Qaintessent.set_creg!(c1, 2)
    Qaintessent.set_creg!(c2, 7)
    @test apply(cgc_ref, ψ) ≈ apply(cgc, ψ)
end

@testset ExtendedTestSet "QRegister" begin
    q1 = qreg(2)
    q2 = qreg(2)
    q3 = qreg(1)
    c1 = creg(2)
    c2 = creg(4)

    cgc = CircuitGateChain([q1,q2,q3])

    N = size(cgc)

    gates = [
        single_qubit_circuit_gate(q1, Y, N),
        single_qubit_circuit_gate(q3[1], Z, N),
        controlled_circuit_gate(q1, (3,), SGate(), N),
        controlled_circuit_gate(q2, 1, TGate(), N),
        controlled_circuit_gate(5, q1, Y, N),
        two_qubit_circuit_gate(q2, q1, SwapGate(), N),
    ]
    cgc(gates)

    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, Y, N),
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(5, Z, N),
        controlled_circuit_gate(1, 3, SGate(), N),
        controlled_circuit_gate(2, 3, SGate(), N),
        controlled_circuit_gate(3, 1, TGate(), N),
        controlled_circuit_gate(4, 1, TGate(), N),
        controlled_circuit_gate(5, 1, Y, N),
        controlled_circuit_gate(5, 2, Y, N),
        two_qubit_circuit_gate(3, 1, SwapGate(), N),
        two_qubit_circuit_gate(4, 2, SwapGate(), N),
    ])

    @test_throws ErrorException("Control and target wires must be disjoint.") CircuitGateChain{N}([
        controlled_circuit_gate(q1, 1, TGate(), N),
    ])

    @test_throws ErrorException("Control and target wires must be disjoint.") CircuitGateChain{N}([
        controlled_circuit_gate(3, q2, TGate(), N),
    ])

    @test_throws ErrorException("Registers used must be of same length") CircuitGateChain{N}([
        controlled_circuit_gate(q1, q3, TGate(), N),
    ])

    ψ = randn(ComplexF64, 2^N)
    @test all(cgc .≈ cgc_ref)
    @test apply(cgc_ref, ψ) ≈ apply(cgc, ψ)
end

@testset ExtendedTestSet "test cgc functor" begin
    c1 = creg(2)
    c2 = creg(4)
    q1 = qreg(3)
    q2 = qreg(3)
    q3 = qreg(4)

    cgc = CircuitGateChain([q1,q2,q3])
    N = size(cgc)

    cgc1 = deepcopy(cgc)
    gates = [
        single_qubit_circuit_gate(2, Y, N),
    ]

    cgc1(gates)
    cgc1_ref = CircuitGateChain{10}(
    [
        single_qubit_circuit_gate(2, Y, N),
    ])
    @test all(cgc1 .≈ cgc1_ref)

    cgc2 = deepcopy(cgc)
    gates = [
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(1, Y, N)
    ]

    cgc2(gates)
    cgc2_ref = CircuitGateChain{10}(
    [
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(1, Y, N)
    ])
    @test all(cgc2.≈ cgc2_ref)

    cgc3 = deepcopy(cgc)
    gates = [
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(q1, Y, N),
        two_qubit_circuit_gate(q1, q2, SwapGate(), N),
        two_qubit_circuit_gate(q1, q3[1], SwapGate(), N),
    ]
    cgc3(gates)

    cgc3_ref = CircuitGateChain{10}(
    [
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(1, Y, N),
        single_qubit_circuit_gate(2, Y, N),
        single_qubit_circuit_gate(3, Y, N),
        two_qubit_circuit_gate(1, 4, SwapGate(), N),
        two_qubit_circuit_gate(2, 5, SwapGate(), N),
        two_qubit_circuit_gate(3, 6, SwapGate(), N),
        two_qubit_circuit_gate(1, 7, SwapGate(), N),
        two_qubit_circuit_gate(2, 7, SwapGate(), N),
        two_qubit_circuit_gate(3, 7, SwapGate(), N),
    ])
    @test all(cgc3 .≈ cgc3_ref)

    cgc4 = deepcopy(cgc)
    gates = [
        controlled_circuit_gate(q1, q2, X, N),
        controlled_circuit_gate(q1, q2[1], Z, N),
        controlled_circuit_gate(q1[3], q2, Y, N),
    ]
    # println(gates)
    cgc4(gates)

    cgc4_ref = CircuitGateChain{10}(
    [
        controlled_circuit_gate(1, 4, X, N),
        controlled_circuit_gate(2, 5, X, N),
        controlled_circuit_gate(3, 6, X, N),
        controlled_circuit_gate(1, 4, Z, N),
        controlled_circuit_gate(2, 4, Z, N),
        controlled_circuit_gate(3, 4, Z, N),
        controlled_circuit_gate(3, 4, Y, N),
        controlled_circuit_gate(3, 5, Y, N),
        controlled_circuit_gate(3, 6, Y, N),
    ])
    @test all(cgc4 .≈ cgc4_ref)

end
