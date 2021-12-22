using Test
using TestSetExtensions
using Qaintessent
using StatsBase
using LinearAlgebra
using CUDA 

if CUDA.functional()
    tol = 1e-6
else
    tol = 1e-15
end

@testset ExtendedTestSet "utility functions" begin
    @testset ExtendedTestSet "util pauli_vector" begin
        n1, n2, n3 = sample(1:15, 3, replace=true)
        @test all(Qaintessent.pauli_vector(n1, n2, n3) .≈ matrix(X) * n1 + matrix(Y) * n2 + matrix(Z) * n3)
    end

    @testset ExtendedTestSet "util intlog2" begin
        @test_throws ErrorException("Logarithm of 0 is undefined") Qaintessent.intlog2(0)
        @test_throws ErrorException("Logarithm of negative number is undefined") Qaintessent.intlog2(-30)
        N = rand(1:50)
        @test Qaintessent.intlog2(2^N) == N
        N = rand(1:200)    
        @test Qaintessent.intlog2(N) <= log(2, N)
        @test Qaintessent.intlog2(N) >= (log(2, N) - 1)
    end

    @testset ExtendedTestSet "util binary_digits" begin
        n1 = abs(rand(Int, 1)[])
        # n1 = 5
        M = Qaintessent.intlog2(n1) + 1
        bd = Qaintessent.binary_digits(M, n1)
        x = 0
        count = 0
        
        for i in bd
            x+= i * 2^count
            count += 1
        end

        @test x ≈ n1
        @test Qaintessent.binary_to_int(bd) ≈ n1
    end

    @testset ExtendedTestSet "util binary_digits!" begin
        # n1 = abs(rand(Int, 1)[])
        n1 = 5
        M = Qaintessent.intlog2(n1) + 1
        m = BitArray{1}(undef, M)
        bd = Qaintessent.binary_digits!(m, n1)
        x = 0
        count = 0
        
        for i in bd
            x+= i * 2^count
            count += 1
        end

        @test x ≈ n1
    end


    @testset ExtendedTestSet "util quaternary_digits" begin
        n1 = abs(rand(Int, 1)[])
        M = Qaintessent.intlog2(n1)
        qd = Qaintessent.quaternary_digits(M, n1)
        x = 0
        count = 0
        for i in qd
            x+= i * 4^count
            count += 1
        end
        @test x ≈ n1
    end

    @testset ExtendedTestSet "util quaternary_digits!" begin
        n1 = abs(rand(Int, 1)[])
        M = Qaintessent.intlog2(n1)
        m = Vector{Int}(undef, M)
        qd = Qaintessent.quaternary_digits!(m, n1)
        x = 0
        count = 0
        for i in qd
            x+= i * 4^count
            count += 1
        end
        @test x ≈ n1
    end

    @testset ExtendedTestSet "util gramm_schmidt!" begin 
        N = rand(4:10)

        @testset "util gramm_schmidt! complex" begin
            a = rand(ComplexQ, (N,N))

            Qaintessent.gramm_schmidt!(a)

            for i in 1:N
                for j in (i+1):N
                    @test norm(dot(a[:, i], a[:, j])) < 2e-6
                end
                @test isapprox(dot(a[:, i], a[:, i]), 1, atol=tol, rtol=tol)
            end
        end

        @testset "util gramm_schmidt! float" begin
            a = rand(FloatQ, (N,N))

            Qaintessent.gramm_schmidt!(a)

            for i in 1:N
                for j in (i+1):N
                    @test norm(dot(a[:, i], a[:, j])) < 2e-6
                end
                @test isapprox(dot(a[:, i], a[:, i]), 1, atol=tol, rtol=tol)
            end
        end
    end

    ##==----------------------------------------------------------------------------------------------------------------------


    @testset ExtendedTestSet "util reduced density matrix" begin
        N = 4
        ψ = randn(ComplexQ, 2^N)
        χ = randn(ComplexQ, 2^N)
        # full density matrix |ψ⟩⟨χ|
        ρ = reshape(kron(conj(χ), ψ), 2^N, 2^N)

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