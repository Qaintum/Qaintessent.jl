using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(cb::CircuitGateChain)
    Qaintessent.matrix(cb) * Qaintessent.matrix(Base.adjoint(cb)) ≈ I
end


@testset ExtendedTestSet "quantum Fourier transform" begin
    for N in 1:4
        c = qft_circuit(N)

        Fmat = [exp(2*π*im*j*k/2^N)/sqrt(2^N) for j in 0:(2^N-1), k in 0:(2^N-1)]
        @test Qaintessent.matrix(c) ≈ Fmat

        @test isunitary(c)
    end
end


@testset ExtendedTestSet "toffoli circuit test" begin
    N = 3
    cgc = toffoli_circuit((1, 2), 3, N)

    # reference
    toffoli = ControlledGate{1,N}(X)

    @test Qaintessent.matrix(toffoli) ≈ Qaintessent.matrix(cgc)
end


@testset ExtendedTestSet "vbe adder test" begin
    for N in 1:4
        M = N * 3 + 1
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
    for N in 1:5
        n = 1
        anc = 0
        while 2^n < N
            anc += N ÷ 2^n - 1
            n += 1
        end
        M = 3N + anc + 1

        cgc = qcla_out_adder_circuit(N)

        a = rand(0:2^(N-1))
        b = rand(0:2^(N-1))

        index = b << N + a
        ψ = fill(0.0 + 0.0 * im, 2^M)

        ψ[index+1] = 1.0

        ψ = apply(cgc, ψ)
        answer = (findall(x -> x == 1, ψ)[1] - 1) >> 2N
        @test answer == a + b
    end
end


@testset ExtendedTestSet "qcla in place adder test" begin
    for N in 1:6
        n = 1
        anc = 0
        while 2^n < N
            anc += N ÷ 2^n - 1
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
        answer = (findall(x -> x == 1, ψ)[1] - 1) >> N
        @test answer == a + b

    end
end
