using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using SparseArrays


##==----------------------------------------------------------------------------------------------------------------------


isunitary(cg::CircuitGate) = (sparse_matrix(cg) * sparse_matrix(Base.adjoint(cg)) ≈ I)


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

    # test sparse_matrix 
    @testset "single qubit circuit gates" begin
        cgs = [circuit_gate(1, X), circuit_gate(2, X), circuit_gate(3, X)]
        m = sparse_matrix(cgs)

        @test m ≈ sparse([8, 7, 6, 5, 4, 3, 2, 1], [1, 2, 3, 4, 5, 6, 7, 8], Float64[1, 1, 1, 1, 1, 1, 1, 1])
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
        @test Qaintessent.sparse_matrix(cg, 3) ≈ kron(Matrix(I, 2, 2), Qaintessent.matrix(Y), Matrix(I, 2, 2))
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
            fill(0, 2, 4) Qaintessent.matrix(HadamardGate()) fill(0, 2, 2);
            fill(0, 2, 6) Qaintessent.matrix(HadamardGate())]
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


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "circuit gate isapprox" begin
    θ = 0.7 * π
    ϕ = 0.4 * π
    n = randn(3); n /= norm(n)
    ϵ = 3*sqrt(eps())
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


##==----------------------------------------------------------------------------------------------------------------------


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
        @test circuit_gate(targ_iwire,    g,  cntrl_iwire)   ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
        @test circuit_gate(targ_iwire,    g, (cntrl_iwire,)) ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
        @test circuit_gate((targ_iwire,), g,  cntrl_iwire)   ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
        @test circuit_gate((targ_iwire,), g, (cntrl_iwire,)) ≈ CircuitGate((targ_iwire, cntrl_iwire), ControlledGate(g, 1))
    end
end
