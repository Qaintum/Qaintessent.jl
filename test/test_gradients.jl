using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


# from https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
function ngradient(f, xs::AbstractArray...)
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in 1:length(x)
        δ = sqrt(eps())
        tmp = x[i]
        x[i] = tmp - δ/2
        y1 = f(xs...)
        x[i] = tmp + δ/2
        y2 = f(xs...)
        x[i] = tmp
        Δ[i] = (y2-y1)/δ
    end
    return grads
end


@testset ExtendedTestSet "circuit gradients" begin

    # construct parametrized circuit
    N = 4
    cgc(θ, ϕ, χ) = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate((1, 4), 2, RzGate(θ), N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(3, PhaseShiftGate(ϕ), N),
        single_qubit_circuit_gate(1, RyGate(χ), N),
    ])
    meas = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N), Hermitian(randn(ComplexF64, 2^N, 2^N))])

    # input quantum state
    ψ = randn(ComplexF64, 2^N)

    # fictitious gradients of cost function with respect to circuit output
    Δ = [0.3, -1.2]

    f(rθ, rϕ, rχ) = dot(Δ, apply(Circuit(cgc(rθ[1], rϕ[1], rχ[1]), meas), ψ))
    θ = 1.5π
    ϕ = 0.3
    χ = √(2)
    @test all(isapprox.(ngradient(f, [θ], [ϕ], [χ]),
        Qaintessent.gradients(Circuit(cgc(θ, ϕ, χ), meas), ψ, Δ), rtol=1e-5, atol=1e-5))
end
