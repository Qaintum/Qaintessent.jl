using Test
using TestSetExtensions
using LinearAlgebra
using Random
using RandomMatrices
using Qaintessent

@testset ExtendedTestSet "unitary compilation" begin
    @testset "compile 1 qubit" begin
        N = 1
        U = Stewart(ComplexF64, 2)
        M = Stewart(ComplexF64, 2)

        cgc = Qaintessent.compile(U, N)

        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled

    end

    @testset "compile diagonal unitaries" begin
        N = 5
        _, U = qr(Stewart(ComplexF64, 2^N))
        M = Stewart(ComplexF64, 2^N)

        cgc = Qaintessent.compile(U, N)

        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "compile 2 qubit unitaries" begin
        N = 2
        U = Stewart(ComplexF64, 2^N)
        M = Stewart(ComplexF64, 2^N)

        cgc = Qaintessent.compile(U, N)
        println(cgc)
        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end
end
