using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "statevector" begin

    @testset "test statevector int constructor" begin
        N = rand(1:10)
        s = Statevector(N)
        @test s.N == N
        @test s.state == zeros(ComplexF64, 2^N)
    end

    @testset "test statevector vector constructor" begin
        N = rand(1:10)
        ψ = randn(ComplexF64, 2^N)
        s = Statevector(ψ)
        @test s.N == N
        @test s.state == ψ
    end

    @testset "test statevector indexing" begin
        N = rand(1:10)
        ψ = randn(ComplexF64, 2^N)
        s = Statevector(ψ)
        r = rand(1:2^N)
        @test s[r] == ψ[r]
    end

    @testset "test statevector attributes" begin
        N = rand(1:10)
        ψ = randn(ComplexF64, 2^N)
        s = Statevector(ψ)
        @test size(s) == size(ψ)
        @test length(s) == length(ψ)
        @test IndexStyle(typeof(s)) == IndexLinear()
    end
end