using Test
using TestSetExtensions
using LinearAlgebra
using Random
using RandomMatrices
using Qaintessent




@testset ExtendedTestSet "compile diagonal unitaries helper functions" begin
    @testset "greyencode" begin
        greyencode.(0:15) == [0 1 3 2 6 7 5 4 12 13 15 14 10 11 9 8]
    end

    @testset "svalue" begin
        greyencode.(0:15) == [0 1 3 2 6 7 5 4 12 13 15 14 10 11 9 8]
    end

    @testset "flip_state" begin
        greyencode.(0:15) == [0 1 3 2 6 7 5 4 12 13 15 14 10 11 9 8]
    end
end

@testset ExtendedTestSet "unitary compilation" begin
    @testset "general compile" begin
        N = 6
        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        M = Stewart(ComplexF64, 2^N)

        cgc = Qaintessent.compile(deepcopy(U), N)

        ψ = rand(ComplexF64, 2^N)
        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "unblocked QR decomposition" begin
        N = 2
        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        Uref = deepcopy(U)
        M = Stewart(ComplexF64, 4)

        QR, τ = Qaintessent.qr_unblocked(deepcopy(U))
        R = diag(QR)
        Q = Matrix{ComplexF64}(I, (2^N, 2^N))
        Id = Matrix{ComplexF64}(I, (2^N, 2^N))

        for i in 1:size(U)[1]-1
            u = pushfirst!(QR[i+1:2^N, i], 1)
            u = u ./ norm(u)
            H = deepcopy(Id)
            H[i:2^N, i:2^N] = H[i:2^N, i:2^N] - 2*u*u'
            Q = H*Q
            U = H*U
        end
        d = diag(U)

        @test diagm(d) ≈ U
        @test diag(U) ≈ R
        @test Q*Uref ≈ U
        @test inv(Q)*U ≈ Uref
    end

    @testset "compile 1 qubit" begin
        N = 1
        U, _ = qr(Stewart(ComplexF64, 2^N))
        U = Matrix(U)
        M = Stewart(ComplexF64, 2)

        cgc = Qaintessent.compile(U, N)

        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "compile diagonal unitaries" begin
        N = 5
        U = diagm(exp.(im .* rand(Float64, 2^N)))
        M = Stewart(ComplexF64, 2^N)

        cgc = Qaintessent.compile(U, N)

        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "compile 2 qubit unitaries" begin
        N = 2
        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        M = Stewart(ComplexF64, 2^N)

        cgc = Qaintessent.compile(U, N)
        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(cgc, ψ)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end
end
