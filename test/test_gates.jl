using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using Random


##==----------------------------------------------------------------------------------------------------------------------


isunitary(g::AbstractGate) = Qaintessent.matrix(g) * Qaintessent.matrix(Base.adjoint(g)) ≈ I


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "quantum gates" begin
    @testset ExtendedTestSet "gate helper functions" begin
        @testset "kron" begin
            @test kron(I,Z) ≈  kron(Matrix(I, (2,2)), matrix(Z))
        end

        @testset "gate Qaintessent.data helper function" begin
            for g in [RxGate, RyGate, RzGate, PhaseShiftGate, EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
                θ = convert(FloatQ, 2π*rand())
                g1 = g(θ)
                @test θ ≈ Qaintessent.data(g1)
            end
            r = rand(FloatQ, 3)
            g = RotationGate(r)
            @test Array(data(g)) ≈ matrix(g)
            @test isnothing(Qaintessent.data(X))
        end

        @testset "gate addition rotation gates" begin
            for g in [RxGate, RyGate, RzGate, PhaseShiftGate, EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
                θ = 2π*rand()
                ϕ = 2π*rand()

                g1 = g(θ)
                g2 = g(ϕ)

                @test g1 + g2 ≈ g(θ+ϕ)
            end
        end

        @testset "gate addition general rotation gate" begin
            
            θ = 2π*rand()
            m = rand(3)
            m = m./norm(m)

            ϕ = 2π*rand()
            n = rand(3)
            n = n./norm(n)

            g1 = RotationGate(θ, m)
            g2 = RotationGate(ϕ, n)

            @test g1 + g2 ≈ RotationGate(m*θ + ϕ*n)
        end

        @testset "gate phase multiplication rotation gates" begin
            for g in [RxGate, RyGate, RzGate, PhaseShiftGate, EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
                θ = 2π*rand()
                α = rand()

                g1 = g(θ)

                @test α*g1 ≈ g(θ*α)
            end
        end

        @testset "gate multiplication general rotation gate" begin
            
            θ = 2π*rand()
            m = rand(3)
            m = m./norm(m)

            α = rand()

            g1 = RotationGate(θ, m)

            @test α*g1 ≈ RotationGate(m*θ*α)
        end
    end

    @testset ExtendedTestSet "gate implementation" begin

        θ = 2π*rand()
        ϕ = 2π*rand()
        n = randn(3); n /= norm(n)

        @testset "basic quantum gates" begin
            for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ), SwapGate(), EntanglementXXGate(θ), EntanglementYYGate(θ), EntanglementZZGate(θ), controlled_not()]
                @test isunitary(g)
                gdag = adjoint(g)

                @test Qaintessent.matrix(gdag) == adjoint(Qaintessent.matrix(g))
                @test LinearAlgebra.ishermitian(g) == (Qaintessent.matrix(gdag) == Qaintessent.matrix(g))
                @test Qaintessent.sparse_matrix(g) == Qaintessent.matrix(g)
            end
        end


        @testset "general rotation gate" begin
            @test Qaintessent.matrix(RotationGate(θ, [1, 0, 0])) ≈ Qaintessent.matrix(RxGate(θ))
            @test Qaintessent.matrix(RotationGate(θ, [0, 1, 0])) ≈ Qaintessent.matrix(RyGate(θ))
            @test Qaintessent.matrix(RotationGate(θ, [0, 0, 1])) ≈ Qaintessent.matrix(RzGate(θ))
            @test Qaintessent.matrix(RotationGate(0, [1, 0, 0])) ≈ ComplexQ[1 0; 0 1]
            @test RotationGate(θ * n) ≈ RotationGate(θ, n)
            @test Qaintessent.matrix(RotationGate(0.3*θ, n)) * Qaintessent.matrix(RotationGate(0.7*θ, n)) ≈ Qaintessent.matrix(RotationGate(θ, n))
            @test Qaintessent.sparse_matrix(RotationGate(θ, n)) ≈ Qaintessent.matrix(RotationGate(θ, n))
        end


        @testset "general rotation gate exceptions" begin
            @test_throws ErrorException("Rotation axis vector must have length 3.") RotationGate(θ, [1, 0, 0, 0])
            @test_throws ErrorException("Norm of rotation axis vector must be 1.")  RotationGate(θ, [1, 2, 0])
            @test_throws ErrorException("Rotation axis vector must have length 3.") RotationGate([0, θ, 0, 1])
        end


        @testset "general unitary gate" begin

            # test MatrixGate
            N = 3
            d = 2
            A = randn(ComplexQ, d^N, d^N)
            U, _ = qr(A)
            U = Array(U)
            gateU = MatrixGate(U)
            gateUdag = adjoint(gateU)
            @test Qaintessent.matrix(gateU) ≈ U
            @test isunitary(gateU)
            @test LinearAlgebra.ishermitian(gateU) == (Qaintessent.matrix(gateUdag) == Qaintessent.matrix(gateU))
            @test LinearAlgebra.ishermitian(gateU) == (Qaintessent.sparse_matrix(gateUdag) == Qaintessent.sparse_matrix(gateU))

            U = Qaintessent.sparse_matrix(X)
            gateU = MatrixGate(U)
            gateUdag = adjoint(gateU)
            @test Qaintessent.sparse_matrix(gateU) ≈ U
            @test isunitary(gateU)
            @test true == (Qaintessent.sparse_matrix(gateUdag) == Qaintessent.sparse_matrix(gateU))
        end

        # not unitary
        @testset "general unitary gate exceptions" begin

            N = 3
            d = 2
            @test_throws ErrorException("Quantum gate must be unitary") MatrixGate(randn(ComplexQ, d^N, d^N))
        end

        @testset "general controlled gate" begin
            # test MatrixGate
            N = 4
            M = 2
            d = 2
            A = randn(ComplexQ, d^M, d^M)
            U, _ = qr(A)
            U = Array(U)
            gateU = MatrixGate(U)
            cgateU = ControlledGate(gateU, M)

            cgateUdag = adjoint(cgateU)
            @test LinearAlgebra.ishermitian(cgateU) == (Qaintessent.sparse_matrix(cgateU) == Qaintessent.sparse_matrix(cgateUdag))
            @test LinearAlgebra.ishermitian(cgateU) == (Qaintessent.matrix(cgateU) == Qaintessent.matrix(cgateUdag))
        end
    end

    @testset ExtendedTestSet "abstract gate tests" begin
        @test_throws Exception adjoint(AbstractGate()) 
        @test_throws Exception num_wires(AbstractGate()) 
        @test_throws Exception matrix(AbstractGate()) 
        @test_throws Exception sparse_matrix(AbstractGate()) 
    end
end