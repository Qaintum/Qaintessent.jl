using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent

@testset ExtendedTestSet "DAG Implementation" begin
    N = 2
    # construct parametrized circuit

    θ = 1.5π
    ϕ = 0.3
    χ = √2
    n = randn(Float64, 3)
    n /= norm(n)
    ωn = 0.2π * n
    M = randn(ComplexF64, 2^N, 2^N)
    M = 0.5*(M + adjoint(M))

    cgc(θ, ϕ, χ, ωn) = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        controlled_circuit_gate((2), 1, Z, N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(2, Z, N),
        single_qubit_circuit_gate(2, Y, N),
    ])
    cgc1 = cgc(θ, ϕ, χ, ωn)
    meas(M) = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N), Hermitian(M)])

    d = Dag(cgc1)

    d = opt_hadamard(d)
    d = opt_adjoint(d)

    cgc = CircuitGateChain(d)

    cgc_refstring =
        "\n" *
        "    1 —[X ]—————————————\n" *
        "        |               \n" *
        "    2 ——•————[Z ]——[Y ]—\n"

    io = IOBuffer()
    show(io, cgc)
    @test String(take!(io)) == cgc_refstring

end

@testset ExtendedTestSet "Optimize VBE Adder" begin
    for N in 1:3
        M = 3N + 1

        # construct parametrized circuit
        adderref = vbe_adder_circuit(N)

        d = Dag(adderref)

        size_new = size(d)
        size_old = Inf

        d = opt_hadamard(d)
        d = opt_adjoint(d)

        adder = CircuitGateChain(d)

        ψ = fill(0.0+0.0*im, 2^M)

        a = abs(rand(Int, 1)[])%(2^N)
        b = abs(rand(Int, 1)[])%(2^N)
        index = b << N + a
        ψ[index+1] = 1.0

        # get reference solution from reference circuit
        ψref = apply(adderref, ψ)
        answerref = ((findall(x->x==1, ψref)[1]-1)%(2^2N)) >> N

        ψsol = apply(adder, ψ)
        answersol = ((findall(x->x==1, ψsol)[1]-1)%(2^2N)) >> N

        @test answersol ≈ answerref
    end
end
