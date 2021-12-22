using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent
using Random

##==----------------------------------------------------------------------------------------------------------------------


"""Checks that tailored apply gives same result as general apply"""

@testset ExtendedTestSet "apply gates to state vector" begin
    N = 6
    θ = 0.7π
    ϕ = 0.4π
    n = randn(3); n /= norm(n)
    ψ = rand(ComplexQ, 2^N)

    @testset "apply basic gates" begin
        # single qubit gates
        for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            i = rand(1:N)
            cg = circuit_gate(i, g)
            cga = CircuitGate{1,AbstractGate}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}

            @test apply(ψ, cg) ≈ apply(ψ, cga)
            @test apply(ψ, cg) ≈ apply(ψ, sparse_matrix(cga, N))
        end
    end

    @testset "apply multiple basic gates" begin
        g = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
        i = randperm(N)[1:3]
        j = randperm(N)[1:3]
        cgs = circuit_gate.(j, g[i])
        cga = CircuitGate[]
        for cg in cgs
            push!(cga, CircuitGate{1,AbstractGate}(cg.iwire, cg.gate))
        end

        @test apply(ψ, cgs) ≈ apply(ψ, cga)
        @test apply(ψ, cgs) ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "apply basic gates exceptions" begin
        cg = CircuitGate[]
        @test_throws ErrorException("Vector of length 0 cannot be applied") apply(ψ, cg)
    end

    @testset "apply moments" begin
        # single qubit gates
        for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            i = rand(1:N)
            cg = Moment([circuit_gate(i, g)])
            cga = Moment([CircuitGate{1,AbstractGate}(cg[1].iwire, cg[1].gate)]) # generate same gate with type AbstractGate{1}

            @test apply(ψ, cg) ≈ apply(ψ, cga)
            @test apply(ψ, cg) ≈ sparse_matrix(cga, N) * ψ
        end
    end

    @testset "apply multiple moments" begin
        g = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
        i = randperm(N)[1:3]
        j = randperm(N)[1:3]
        cgs = circuit_gate.(j, g[i])
        cga = CircuitGate[]
        for cg in cgs
            push!(cga, CircuitGate{1,AbstractGate}(cg.iwire, cg.gate))
        end

        @test apply(ψ, Moment(cgs)) ≈ apply(ψ, Moment(cga))
        @test apply(ψ, Moment(cgs)) ≈ sparse_matrix(cga, N) * ψ
    end


    @testset "apply moments exceptions" begin
        m = Moment[]
        @test_throws ErrorException("Vector of length 0 cannot be applied") apply(ψ, m)
    end


    @testset "apply basic controlled gates" begin
        # control gate
        for g in [X, Y, Z, RotationGate(θ,n), PhaseShiftGate(ϕ)]
            i = rand(1:N)
            j = rand([1:i-1; i+1:N])
            cg = circuit_gate(i, g, j)
            cga = CircuitGate{2,ControlledGate{typeof(g)}}(cg.iwire, cg.gate)

            @test apply(ψ, cg) ≈ apply(ψ, cga)
            @test apply(ψ, cg) ≈ sparse_matrix(cga, N) * ψ
        end
    end

    @testset "apply swap gate" begin
        i = rand(1:N)
        j = rand([1:i-1; i+1:N])
        cg = circuit_gate((i, j), SwapGate())
        cga = CircuitGate{2,SwapGate}(cg.iwire, cg.gate)

        @test apply(ψ, cg) ≈ apply(ψ, cga)
        @test apply(ψ, cg) ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "apply 1-qubit MatrixGate" begin
        # MatrixGate: one qubit
        d = 2
        A = rand(ComplexQ, d, d)
        U, R = qr(A)
        U = Array(U);
        g = MatrixGate(U)
        i = rand(1:N)
        cg = circuit_gate(i, g)
        cga = CircuitGate{1,MatrixGate}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}

        @test apply(ψ, cg) ≈ apply(ψ, cga)
        @test apply(ψ, cg) ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "apply k-qubit MatrixGate" begin
        # MatrixGate: k qubits
        k = rand(1:N)
        A = rand(ComplexQ, 2^k, 2^k)
        U, R = qr(A)
        U = Array(U);
        g = MatrixGate(U)
        iwire = [rand(1:N)]
        for j in 1:k-1
            l = rand(1:N-j)
            i = setdiff([1:N...], iwire)[l]
            push!(iwire, i)
        end
        sort!(iwire)
        cga = CircuitGate{k,MatrixGate}((iwire...,), g)
        m = sparse_matrix(cga, N)
        @test apply(ψ, cga) ≈ m*ψ
    end
end