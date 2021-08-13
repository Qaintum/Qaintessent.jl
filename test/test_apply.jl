using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


"""Checks that tailored apply gives same result as general apply"""

@testset ExtendedTestSet "apply gates to state vector" begin
    N = 6
    θ = 0.7π
    ϕ = 0.4π
    n = randn(3); n /= norm(n)
    ψ = rand(ComplexF64, 2^N)

    @testset "apply basic gates" begin
        # single qubit gates
        for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            i = rand(1:N)
            cg = circuit_gate(i, g)
            cga = CircuitGate{1,AbstractGate}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}

            @test apply(ψ, cg) ≈ apply(ψ, cga)
            @test apply(ψ, cg) ≈ sparse_matrix(cga, N) * ψ
        end
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
        A = rand(ComplexF64, d, d)
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
        A = rand(ComplexF64, 2^k, 2^k)
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
