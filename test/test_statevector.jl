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
        @test s.state == zeros(ComplexQ, 2^N)
    end

    @testset "test statevector vector constructor" begin
        N = rand(1:10)
        ψ = randn(ComplexQ, 2^N)
        s = Statevector(ψ)
        @test s.N == N
        @test s.state == ψ

        if ComplexQ == ComplexF64
            ComplexTest = ComplexF32
        else
            ComplexTest = ComplexF64
        end

        ψ = randn(ComplexTest, 2^N)
        s = Statevector(ψ)
        @test typeof(s.state) == Vector{ComplexQ}
    end

    @testset "test statevector indexing" begin
        N = rand(1:10)
        ψ = randn(ComplexQ, 2^N)
        s = Statevector(ψ)
        r = rand(1:2^N)
        @test s[r] == ψ[r]
    end

    @testset "test statevector attributes" begin
        N = rand(1:10)
        ψ = randn(ComplexQ, 2^N)
        s = Statevector(ψ)
        @test size(s) == size(ψ)
        @test length(s) == length(ψ)
        @test IndexStyle(typeof(s)) == IndexLinear()
    end

    @testset ExtendedTestSet "util reduced density matrix for statevector" begin
        N = 4
        ψ = Statevector(randn(ComplexQ, 2^N))
        χ = Statevector(randn(ComplexQ, 2^N))
        # full density matrix |ψ⟩⟨χ|
        ρ = reshape(kron(conj(χ.state), ψ.state), 2^N, 2^N)
        id = Matrix{ComplexQ}(I, 2, 2)
        A = randn(ComplexQ, 2, 2)
        B = randn(ComplexQ, 2, 2)
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

    @testset ExtendedTestSet "util reduced density matrix for statevector" begin
        N = 4
        ψ = Statevector(randn(ComplexQ, 2^N))
        χ = Statevector(randn(ComplexQ, 2^N))
        # full density matrix |ψ⟩⟨χ|
        ρ = reshape(kron(conj(χ.state), ψ.state), 2^N, 2^N)
        id = Matrix{ComplexQ}(I, 2, 2)
        A = randn(ComplexQ, 2, 2)
        B = randn(ComplexQ, 2, 2)
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


end