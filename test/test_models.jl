using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(cb::Vector{CircuitGate})
    matrix(cb) * matrix(Base.adjoint(cb)) ≈ I
end


@testset ExtendedTestSet "quantum Fourier transform" begin
    for N in 1:4
        c = qft_circuit(N)

        Fmat = [exp(2*π*im*j*k/2^N)/sqrt(2^N) for k in 0:(2^N-1), j in 0:(2^N-1)]

        @test Qaintessent.matrix(c) ≈ Fmat
        @test isunitary(c)
    end
end


@testset ExtendedTestSet "toffoli circuit test" begin
    N = 3
    cgc = toffoli_circuit(1, (3, 2))

    # reference
    toffoli = circuit_gate(1, X, 3, 2)

    @test Qaintessent.matrix(toffoli) ≈ Qaintessent.matrix(cgc)
end


@testset ExtendedTestSet "vbe adder test" begin
    for N in 1:4
        M = 3N + 1
        adder = vbe_adder_circuit(N)
        ψ = fill(0.0 + 0.0im, 2^M)
        a = rand(0:2^(N-1))
        b = rand(0:2^(N-1))

        index = b << N + a
        ψ[index+1] = 1.0
        ψ = apply(adder, ψ)
        answer = ((findall(x -> x == 1, ψ)[1] - 1) % (2^2N)) >> N
        @test answer == (a + b) % (2^N)
    end
end


@testset ExtendedTestSet "qcla out of place adder test" begin
    for N in 1:4
        anc = N - count_ones(N) - floor(Int, log(N))
        M = 3N + anc + 1

        cgc = qcla_out_adder_circuit(N)
        ψ = fill(0.0 + 0.0 * im, 2^M)
        a = rand(0:2^(N-1))
        b = rand(0:2^(N-1))

        index = b << N + a

        ψ[index+1] = 1.0

        ψ = apply(cgc, ψ)
        answer = (findall(x -> x == 1, ψ)[1] - 1) >> 2N
        @test answer == a + b
    end
end


@testset ExtendedTestSet "qcla in place adder test" begin
    for N in 1:4
        n = 1
        anc = 0
        while 2^n < N
            anc += N ÷ 2^n -1
            n += 1
        end
        M = 3N + anc
        cgc = qcla_inplace_adder_circuit(N)
        a = rand(0:2^(N-1))
        b = rand(0:2^(N-1))

        index = b << N + a
        ψ = fill(0.0 + 0.0 * im, 2^M)

        ψ[index+1] = 1.0

        ψ = apply(cgc, ψ)
        answer = ((findall(x -> x == 1, ψ)[1] - 1) >> N)
        @test answer == a + b
    end
end
