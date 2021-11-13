using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


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

    @testset "moments approx" begin
        g = CircuitGate((2,), ZGate())
        h = CircuitGate((1,), SGate())
        m = Moment([g, h])
        m1 = Moment([g, h])
        m2 = Moment([h, g])

        @test !(m1 ≈ m2)
        @test m1 ≈ m
    end

    @testset "moments array operations" begin
        m = Moment([CircuitGate((2,), ZGate()), CircuitGate((1,), SGate())])

        @test all(m[:] .≈ m.gates)
        @test all(m[1:2] .≈ m.gates)
        @test lastindex(m) == length(m.gates)
        @test isempty(m) == false
        @test CircuitGate((1,), SGate()) ≈ pop!(m)
        @test isempty(m) == false
        @test CircuitGate((2,), ZGate()) ≈ pop!(m)
        @test isempty(m) == true
        @test firstindex(m) == 1
        @test_throws ArgumentError("array must be non-empty") pop!(m)
    end

    @testset "moments sparse_matrix" begin
        g = CircuitGate((2,), ZGate())
        h = CircuitGate((1,), SGate())
        m = Moment([g, h])
        ms = Moment[Moment([g]), Moment([h])]
        
        @test sparse_matrix(m) ≈ sparse_matrix(g, 2) * sparse_matrix(h, 2)
        @test sparse_matrix(ms) ≈ sparse_matrix(g, 2) * sparse_matrix(h, 2)
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

        @test all(reverse!(m) .≈ m_reverse)
        @test m[1] ≈ h
        @test m[2] ≈ g
    end

    @testset "moments exceptions" begin
        @testset "sparse_matrix for empty Moments" begin
            @test_throws ErrorException("Vector of length 0 cannot be converted to matrix") sparse_matrix(Moment[])
        end
    end
end