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