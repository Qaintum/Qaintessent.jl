using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


# adapted from https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
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
        if eltype(x) <: Complex
            # derivative with respect to imaginary part
            x[i] = tmp - im*δ/2
            y1 = f(xs...)
            x[i] = tmp + im*δ/2
            y2 = f(xs...)
            x[i] = tmp
            Δ[i] += im*(y2-y1)/δ
        end
    end
    return grads
end


@testset ExtendedTestSet "circuit gradients" begin

    # construct parametrized circuit
    N = 4
    cgc(θ, ϕ, χ, ωn) = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate((1, 4), 2, RzGate(θ), N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(3, PhaseShiftGate(ϕ), N),
        single_qubit_circuit_gate(3, RotationGate(ωn), N),
        single_qubit_circuit_gate(1, RyGate(χ), N),
    ])
    # measurement operators
    meas(M) = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N), Hermitian(M)])

    # parameter values
    θ = 1.5π
    ϕ = 0.3
    χ = √2
    n = randn(Float64, 3)
    n /= norm(n)
    ωn = 0.2π * n
    M = randn(ComplexF64, 2^N, 2^N)
    M = 0.5*(M + adjoint(M))

    c = Circuit(cgc(θ, ϕ, χ, ωn), meas(M))

    # input quantum state
    ψ = randn(ComplexF64, 2^N)

    # fictitious gradients of cost function with respect to circuit output
    Δ = [0.3, -1.2]

    dc = Qaintessent.gradients(c, ψ, Δ)[1]

    f(rθ, rϕ, rχ, ωn, M) = dot(Δ, apply(Circuit(cgc(rθ[], rϕ[], rχ[], ωn), meas(M)), ψ))
    # numeric gradients
    ngrad = ngradient(f, [θ], [ϕ], [χ], ωn, M)
    # symmetrize gradient with respect to measurement operator
    ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
    @test all(isapprox.(ngrad,
        (dc.cgc[2][1].gate.U.θ,
         dc.cgc[4][1].gate.ϕ,
         dc.cgc[6][1].gate.θ,
         dc.cgc[5][1].gate.nθ,
         dc.meas.mops[2]), rtol=1e-5, atol=1e-5))
end


@testset ExtendedTestSet "circuit gradients" begin

    # construct parametrized circuit
    N = 4
    cgc(θ, ϕ, χ, ωn) = CircuitGateChain{N}([
        Moment{N}(
        [single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate((1, 4), 2, RzGate(θ), N)]
        ),
        Moment{N}(two_qubit_circuit_gate(2, 3, SwapGate(), N)),
        Moment{N}(single_qubit_circuit_gate(3, PhaseShiftGate(ϕ), N)),
        Moment{N}(
        [single_qubit_circuit_gate(3, RotationGate(ωn), N),
        single_qubit_circuit_gate(1, RyGate(χ), N)]
        ),
    ])
    # measurement operators
    meas(M) = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N), Hermitian(M)])

    # parameter values
    θ = 1.5π
    ϕ = 0.3
    χ = √2
    n = randn(Float64, 3)
    n /= norm(n)
    ωn = 0.2π * n
    M = randn(ComplexF64, 2^N, 2^N)
    M = 0.5*(M + adjoint(M))

    c = Circuit(cgc(θ, ϕ, χ, ωn), meas(M))

    # input quantum state
    ψ = randn(ComplexF64, 2^N)

    # fictitious gradients of cost function with respect to circuit output
    Δ = [0.3, -1.2]

    dc = Qaintessent.gradients(c, ψ, Δ)[1]

    f(rθ, rϕ, rχ, ωn, M) = dot(Δ, apply(Circuit(cgc(rθ[], rϕ[], rχ[], ωn), meas(M)), ψ))
    # numeric gradients
    ngrad = ngradient(f, [θ], [ϕ], [χ], ωn, M)
    # symmetrize gradient with respect to measurement operator
    ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
    @test all(isapprox.(ngrad,
        (dc.cgc[1][2].gate.U.θ,
         dc.cgc[3][1].gate.ϕ,
         dc.cgc[4][2].gate.θ,
         dc.cgc[4][1].gate.nθ,
         dc.meas.mops[2]), rtol=1e-5, atol=1e-5))
end
