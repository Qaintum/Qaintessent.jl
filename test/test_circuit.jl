using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(cg::CircuitGate)
    sparse_matrix(cg) * sparse_matrix(Base.adjoint(cg)) ≈ I
end

@testset ExtendedTestSet "circuit gates" begin
    θ = 0.7 * π
    ϕ = 0.4 * π
    n = randn(3); n /= norm(n)

    # single qubit gates
    @testset "single qubit circuit gates" begin
        for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            cg = CircuitGate((2,), g)
            cgadj = adjoint(cg)
            Qaintessent.sparse_matrix(cgadj.gate) == adjoint(Qaintessent.sparse_matrix(cg.gate))
            @test LinearAlgebra.ishermitian(cg) == (Qaintessent.sparse_matrix(cg) == Qaintessent.sparse_matrix(adjoint(cg)))
        end
    end

    # two qubit gates
    @testset "two qubit circuit gates" begin
        for g in [EntanglementXXGate(θ), EntanglementYYGate(θ), EntanglementZZGate(θ), controlled_not(), SwapGate()]
            cg = CircuitGate((2, 3), g)
            cgadj = adjoint(cg)
            Qaintessent.sparse_matrix(cgadj.gate) == adjoint(Qaintessent.sparse_matrix(cg.gate))
            @test LinearAlgebra.ishermitian(cg) == (Qaintessent.sparse_matrix(cg) == Qaintessent.sparse_matrix(adjoint(cg)))
        end
    end

    # Y acting on second wire
    @testset "apply circuit gate to second wire" begin
        cg = CircuitGate((2,), Y)
        @test Qaintessent.sparse_matrix(cg, 3) ≈ kron(kron(Matrix(I, 2, 2), Qaintessent.sparse_matrix(Y)), Matrix(I, 2, 2))
        @test isunitary(cg)
    end

    # flip control and target
    @testset "flip control and target circuit gate" begin
        cg = CircuitGate((2, 1), controlled_not())
        @test Qaintessent.sparse_matrix(cg) ≈ [1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]
        @test isunitary(cg)
    end

    # third qubit as control and first qubit as target
    @testset "shift control and target circuit gate" begin
        cg = circuit_gate(1, HadamardGate(), 3)
        @test Qaintessent.sparse_matrix(cg) ≈ [
            Matrix(I, 4, 4) fill(0, 4, 2) fill(0, 4, 2);
            fill(0, 2, 4) Qaintessent.sparse_matrix(HadamardGate()) fill(0, 2, 2);
            fill(0, 2, 6) Qaintessent.sparse_matrix(HadamardGate())]
        @test isunitary(cg)
    end

    @testset "circuit gate exceptions" begin
        H = HadamardGate()
        S = SwapGate()

        N = 2
        @test_throws ErrorException("SwapGate affects 2 wires but 0 wires, (), were passed.") CircuitGate{0,SwapGate}(NTuple{0,Int}(), S)
        @test_throws ErrorException("Wire indices must be unique.") CircuitGate{2,SwapGate}((1, 1), S)
        @test_throws ErrorException("Wire index cannot be smaller than 1.") CircuitGate{2,SwapGate}((1, -1), S)
    end
end

@testset ExtendedTestSet "circuit gate isapprox" begin
    θ = 0.7 * π
    ϕ = 0.4 * π
    n = randn(3); n /= norm(n)
    ϵ = 3sqrt(eps())
    sqg = [RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
    sqḡ = [RxGate(θ + eps()), RyGate(θ + eps()), RzGate(θ + eps()), RotationGate(θ + eps(), n), PhaseShiftGate(ϕ + eps())]
    sqĝ = [RxGate(θ + ϵ), RyGate(θ + ϵ), RzGate(θ + ϵ), RotationGate(θ + ϵ, n), PhaseShiftGate(ϕ + ϵ)]

    for (i, g) in enumerate(sqg)
        cg1 = CircuitGate((2,), sqg[i])
        cg2 = CircuitGate((2,), sqḡ[i])
        cg3 = CircuitGate((2,), sqĝ[i])
        
        @test cg1 ≈ cg2
        @test !(cg1 ≈ cg3)
    end
end

@testset ExtendedTestSet "circuit gate helper functions" begin
    N = 5
    @testset ExtendedTestSet "circuit gate single qubit helper function" begin
        iwire = rand(1:N)
        g = XGate()
        @test circuit_gate(iwire, g) ≈ CircuitGate((iwire,), g)
    end

    @testset "circuit gate two qubit helper function" begin
        iwire1 = rand(1:N)
        iwire2 = rand(vcat(1:iwire1 - 1..., iwire1 + 1:N...))
        g2 = SwapGate()
        @test circuit_gate(iwire1, iwire2, g2) ≈ CircuitGate((iwire1, iwire2), g2)
    end

    @testset "circuit gate controlled gate helper function" begin
        cntrl_iwire = rand(1:N)
        targ_iwire = rand(vcat(1:cntrl_iwire - 1..., cntrl_iwire + 1:N...))
        g = YGate()
        @test circuit_gate(targ_iwire,     g,   cntrl_iwire) ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
        @test circuit_gate(targ_iwire,     g, (cntrl_iwire,)) ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
        @test circuit_gate((targ_iwire,),   g, cntrl_iwire) ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
        @test circuit_gate((targ_iwire,),   g, (cntrl_iwire,)) ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
    end
end

@testset ExtendedTestSet "moments" begin
    N = 3
    @testset "moments constructor" begin
        g = CircuitGate((2,), XGate())
        m = Moment(g)
        @test m[1] === g

        h = CircuitGate((3,), YGate())
        m = Moment([g,h])
        @test m[1] === g
        @test m[2] === h
    end

    @testset "moments constructor exceptions" begin
        g = CircuitGate((2,), XGate())
        h = CircuitGate((2,), YGate())
        @test_throws ErrorException("Only gates on different wires are allowed in a Moment") cgc_ref =
            Moment([g, h])
    end

    @testset "moments adjoint" begin
        g = CircuitGate((2,), ZGate())
        h = CircuitGate((1,), SGate())
        m = Moment([g, h])
        madj = adjoint(m)
        @test madj[1] ≈ adjoint(h)
        @test madj[2] ≈ adjoint(g)

        @test sparse_matrix(adjoint(m)) ≈ adjoint(sparse_matrix(m))
    end

    @testset "moments reverse" begin
        g = CircuitGate((2,), ZGate())
        h = CircuitGate((1,), SGate())
        m = Moment([g, h])
        m_reverse = reverse(m)
        @test m_reverse[2] ≈ g
        @test m_reverse[1] ≈ h
        @test m[1] ≈ g
        @test m[2] ≈ h
    end
end

@testset ExtendedTestSet "measurement operator construction" begin
    N = 5
    meas1 = MeasurementOperator.([3 * Matrix{Float64}(I, 2^N, 2^N), 2 * Matrix{Float64}(I, 2^N, 2^N)], ((1,2,3,4,5),))
    @test Qaintessent.check_commute(meas1) == true

    meas2 = MeasurementOperator.([X, Y], [(1,),(2,)])
    @test Qaintessent.check_commute(meas2) == true
end

@testset ExtendedTestSet "measurement operator abstract matrix exceptions" begin
    N = 2
    small_matrix = Matrix{Float64}(I, 2^(N)-1, 2^(N)-1)
    large_matrix = Matrix{Float64}(I, 2^(N)+1, 2^(N)+1)
    uneven_matrix = Matrix{Float64}(I, 2^(N), 2^(N+1))
    iwires = (3,4)
    non_hermitian_matrix = Qaintessent.sparse_matrix(circuit_gate(1, RxGate(0.2), 2))

    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(small_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(small_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(large_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(large_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(uneven_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(uneven_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be Hermitian.") MeasurementOperator(non_hermitian_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be Hermitian.") MeasurementOperator(non_hermitian_matrix, iwires)

    noncomm_matrix_1 = sparse_matrix(circuit_gate(1, HadamardGate(), 2))
    noncomm_matrix_2 = sparse_matrix(circuit_gate(1, X, 2))
    m = MeasurementOperator.([noncomm_matrix_1, noncomm_matrix_2], (iwires,))
    @test Qaintessent.check_commute(m) == false
end

@testset ExtendedTestSet "measurement operator abstract matrix exceptions" begin
    iwires = (1,)

    non_hermitian_matrix = RxGate(0.2)
    @test_throws ErrorException("Measurement operator must be Hermitian.") MeasurementOperator(non_hermitian_matrix, iwires)
    iwires = (1,)
    noncomm_gate_1 = HadamardGate()
    noncomm_gate_2 = X
    m = MeasurementOperator.([noncomm_gate_1, noncomm_gate_2], (iwires,))
    @test Qaintessent.check_commute(m) == false
end

@testset ExtendedTestSet "circuit construction" begin
    N = 1
    cgs = [
        circuit_gate(1, X),
        circuit_gate(1, HadamardGate()),
        circuit_gate(1, Z),
        circuit_gate(1, Y),
        ]
    ref_moments = Moment.(cgs)
    meas = MeasurementOperator(Matrix{Float64}(I, 2^N, 2^N), (1,))
    c = Circuit{N}(cgs, [meas])
    for i in 1:length(cgs)
        @test c[i] ≈ ref_moments[i]
    end
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.moments, ψ)

    ψs = apply(cgs, ψ)
    ψref = apply(c, ψ)
    
    @test [dot(ψs, apply(meas,ψs))] ≈ apply(c, ψ)
end

@testset ExtendedTestSet "circuit construction with circuit gate measurements" begin
    N = 3
    cgc = [
        circuit_gate(1, X),
        circuit_gate(2, HadamardGate()),
        circuit_gate(3, Z),
        circuit_gate(1, Y),
    ]
    iwires = [(1,), (2,), (3,)]
    meas = MeasurementOperator.([X, X, X], iwires)
    c = Circuit{N}(cgc, meas)
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.moments, ψ)

    ψs = apply(cgc, ψ)
    @test [dot(ψs, apply(m, ψs)) for m in meas] ≈ apply(c, ψ)
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
    @testset ExtendedTestSet "reduced density matrix correctness" begin
        @test sum(kron(A, B) .* rdm(N, (4, 2), ψ, χ)) ≈ sum(kron(B, id, A, id) .* ρ)
    end

    @testset ExtendedTestSet "reduced density matrix exceptions" begin
        @test_throws ErrorException("Need at least one wire to act on.") rdm(N, (), ψ, χ)
        @test_throws ErrorException("Number of gate wires cannot be larger than total number of wires.") rdm(N, (1, 2, 3, 4, 5), ψ, χ)
        @test_throws ErrorException("Wire indices must be unique.") rdm(N, (2, 2), ψ, χ)
        @test_throws ErrorException("Wire index cannot be smaller than 1.") rdm(N, (-1, 2), ψ, χ)
        @test_throws ErrorException("Wire index cannot be larger than total number of wires.") rdm(N, (5, 1), ψ, χ)
    end
end


# @testset ExtendedTestSet "CRegister" begin
#     N = 5

#     c1 = creg(2)
#     c2 = creg(4)

#     cgc = CircuitGateChain{N}([c1, c2])

#     Qaintessent.set_creg!(c1, 2)
#     Qaintessent.set_creg!(c2, 1, true)

#     gates = [
#         single_qubit_circuit_gate(2, Y, N),
#         controlled_circuit_gate((2,), reg_check(c2, 3), SGate(), N),
#         controlled_circuit_gate(2, (c2[2], 3), YGate(), N),
#         two_qubit_circuit_gate(2, 3, SwapGate(), N),
#         single_qubit_circuit_gate(5, RxGate(1.5π), N),
#         single_qubit_circuit_gate(3, RyGate(1.5π), N),
#         controlled_circuit_gate(2, c2[1], ZGate(), N),
#     ]
#     append!(cgc, gates)

#     cgc_ref = CircuitGateChain{N}([
#         single_qubit_circuit_gate(2, Y, N),
#         two_qubit_circuit_gate(2, 3, SwapGate(), N),
#         single_qubit_circuit_gate(5, RxGate(1.5π), N),
#         single_qubit_circuit_gate(3, RyGate(1.5π), N),
#         single_qubit_circuit_gate(2, ZGate(), N)
#     ])
#     ψ = randn(ComplexF64, 2^N)
#     @test apply(cgc_ref, ψ) ≈ apply(cgc, ψ)

#     c3 = qreg(4)

# end

# @testset ExtendedTestSet "test add_creg!" begin
#     N = 7
#     c1 = creg(2)
#     c2 = creg(4)
#     c3 = creg(3)

#     cgc = CircuitGateChain{N}()

#     add_creg!(cgc, c1)
#     add_creg!(cgc, c2)

#     @test_throws ErrorException("CRegister is already used in another CircuitGateChain object") add_creg!(cgc, c2)
#     @test_throws ErrorException("Register object has yet to be used in a CircuitGateChain object") gates = [controlled_circuit_gate(c3[1], 2, SGate(), N),]

#     gates = [
#         single_qubit_circuit_gate(2, Y, N),
#         controlled_circuit_gate(2, reg_check(c2, 7), SGate(), N),
#         controlled_circuit_gate(2, (reg_check(c2, 3), 3), YGate(), N),
#         two_qubit_circuit_gate(2, 3, SwapGate(), N),
#         single_qubit_circuit_gate(5, RxGate(1.5π), N),
#         single_qubit_circuit_gate(3, RyGate(1.5π), N),
#         controlled_circuit_gate(2, (c2[1],), ZGate(), N),
#     ]
#     append!(cgc, gates)

#     cgc_ref = CircuitGateChain{N}([
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(2, SGate(), N),
#         two_qubit_circuit_gate(2, 3, SwapGate(), N),
#         single_qubit_circuit_gate(5, RxGate(1.5π), N),
#         single_qubit_circuit_gate(3, RyGate(1.5π), N),
#         single_qubit_circuit_gate(2, ZGate(), N)
#     ])

#     ψ = randn(ComplexF64, 2^N)
#     @test !(apply(cgc_ref, ψ) ≈ apply(cgc, ψ))
#     Qaintessent.set_creg!(c1, 2)
#     Qaintessent.set_creg!(c2, 7)
#     @test apply(cgc_ref, ψ) ≈ apply(cgc, ψ)
# end

# @testset ExtendedTestSet "QRegister" begin
#     q1 = qreg(2)
#     q2 = qreg(2)
#     q3 = qreg(1)
#     c1 = creg(2)
#     c2 = creg(4)

#     cgc = CircuitGateChain([q1,q2,q3])

#     N = size(cgc)

#     gates = [
#         single_qubit_circuit_gate(q1, Y, N),
#         single_qubit_circuit_gate(q3[1], Z, N),
#         controlled_circuit_gate((3,), q1, SGate(), N),
#         controlled_circuit_gate(1, q2, TGate(), N),
#         controlled_circuit_gate(q1, 5, Y, N),
#         two_qubit_circuit_gate(q2, q1, SwapGate(), N),
#     ]
#     append!(cgc, gates)

#     cgc_ref = CircuitGateChain{N}([
#         single_qubit_circuit_gate(1, Y, N),
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(5, Z, N),
#         controlled_circuit_gate(3, 1, SGate(), N),
#         controlled_circuit_gate(3, 2, SGate(), N),
#         controlled_circuit_gate(1, 3, TGate(), N),
#         controlled_circuit_gate(1, 4, TGate(), N),
#         controlled_circuit_gate(1, 5, Y, N),
#         controlled_circuit_gate(2, 5, Y, N),
#         two_qubit_circuit_gate(3, 1, SwapGate(), N),
#         two_qubit_circuit_gate(4, 2, SwapGate(), N),
#     ])

#     @test_throws ErrorException("Control and target wires must be disjoint.") CircuitGateChain{N}([
#         controlled_circuit_gate(q1, 1, TGate(), N),
#     ])

#     @test_throws ErrorException("Control and target wires must be disjoint.") CircuitGateChain{N}([
#         controlled_circuit_gate(3, q2, TGate(), N),
#     ])

#     @test_throws ErrorException("Registers used must be of same length") CircuitGateChain{N}([
#         controlled_circuit_gate(q1, q3, TGate(), N),
#     ])

#     ψ = randn(ComplexF64, 2^N)
#     @test all(cgc .≈ cgc_ref)
#     @test apply(cgc_ref, ψ) ≈ apply(cgc, ψ)
# end

# @testset ExtendedTestSet "cgc appending gates" begin
#     c1 = creg(2)
#     c2 = creg(4)
#     q1 = qreg(3)
#     q2 = qreg(3)
#     q3 = qreg(4)

#     cgc = CircuitGateChain([q1,q2,q3])
#     N = size(cgc)

#     cgc1 = deepcopy(cgc)
#     gates = [
#         single_qubit_circuit_gate(2, Y, N),
#     ]

#     append!(cgc1, gates)
#     cgc1_ref = CircuitGateChain{10}(
#     [
#         single_qubit_circuit_gate(2, Y, N),
#     ])
#     @test all(cgc1 .≈ cgc1_ref)

#     cgc2 = deepcopy(cgc)
#     gates = [
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(1, Y, N)
#     ]

#     append!(cgc2, gates)
#     cgc2_ref = CircuitGateChain{10}(
#     [
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(1, Y, N)
#     ])
#     @test all(cgc2.≈ cgc2_ref)

#     cgc3 = deepcopy(cgc)
#     gates = [
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(q1, Y, N),
#         two_qubit_circuit_gate(q1, q2, SwapGate(), N),
#         two_qubit_circuit_gate(q1, q3[1], SwapGate(), N),
#     ]
#     append!(cgc3, gates)

#     cgc3_ref = CircuitGateChain{10}(
#     [
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(1, Y, N),
#         single_qubit_circuit_gate(2, Y, N),
#         single_qubit_circuit_gate(3, Y, N),
#         two_qubit_circuit_gate(1, 4, SwapGate(), N),
#         two_qubit_circuit_gate(2, 5, SwapGate(), N),
#         two_qubit_circuit_gate(3, 6, SwapGate(), N),
#         two_qubit_circuit_gate(1, 7, SwapGate(), N),
#         two_qubit_circuit_gate(2, 7, SwapGate(), N),
#         two_qubit_circuit_gate(3, 7, SwapGate(), N),
#     ])
#     @test all(cgc3 .≈ cgc3_ref)

#     cgc4 = deepcopy(cgc)
#     gates = [
#         controlled_circuit_gate(q1, q2, X, N),
#         controlled_circuit_gate(q1, q2[1], Z, N),
#         controlled_circuit_gate(q1[3], q2, Y, N),
#     ]
#     # println(gates)
#     append!(cgc4, gates)

#     cgc4_ref = CircuitGateChain{10}(
#     [
#         controlled_circuit_gate(1, 4, X, N),
#         controlled_circuit_gate(2, 5, X, N),
#         controlled_circuit_gate(3, 6, X, N),
#         controlled_circuit_gate(1, 4, Z, N),
#         controlled_circuit_gate(2, 4, Z, N),
#         controlled_circuit_gate(3, 4, Z, N),
#         controlled_circuit_gate(3, 4, Y, N),
#         controlled_circuit_gate(3, 5, Y, N),
#         controlled_circuit_gate(3, 6, Y, N),
#     ])
#     @test all(cgc4 .≈ cgc4_ref)
# end
