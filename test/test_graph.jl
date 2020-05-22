using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent

@testset ExtendedTestSet "circuit gradients" begin
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
        two_qubit_circuit_gate(1, 2, IIGate(), N),
        two_qubit_circuit_gate(1, 2, IIGate(), N),
        controlled_circuit_gate((2), 1, Z, N),
        two_qubit_circuit_gate(1, 2, IIGate(), N),
        two_qubit_circuit_gate(1, 2, IIGate(), N),
        single_qubit_circuit_gate(1, HadamardGate(), N),
        single_qubit_circuit_gate(1, X, N),
        single_qubit_circuit_gate(2, Z, N),
        single_qubit_circuit_gate(2, Y, N),
    ])
    cgc1 = cgc(θ, ϕ, χ, ωn)
    println(cgc1)
    meas(M) = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N), Hermitian(M)])

    d = Dag(cgc1)

    println(cgc2)
    d = opt_hadamard(d)
    d = opt_adjoint(d)
    println("")
    println(d)

    cgc = CircuitGateChain(d)
    println(cgc)

end
