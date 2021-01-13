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


@testset ExtendedTestSet "gate gradients density matrix" begin

    @testset "single qubit gates" begin

        # fictitious density matrix
        ρ = DensityMatrix(randn(Float64, 4), 1)

        # fictitious gradients of cost function with respect to output density matrix
        Δ = DensityMatrix(randn(Float64, 4), 1)

        for g in [RxGate, RyGate, RzGate]
            f(θ) = dot(Δ.v, apply(CircuitGate((1,), g(θ[])), ρ).v)
            θ = 2π*rand()
            ngrad = ngradient(f, [θ])
            dg = Qaintessent.backward_density(g(θ), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)
        end

        begin
            f(ϕ) = dot(Δ.v, apply(CircuitGate((1,), PhaseShiftGate(ϕ[])), ρ).v)
            ϕ = 2π*rand()
            ngrad = ngradient(f, [ϕ])
            dg = Qaintessent.backward_density(PhaseShiftGate(ϕ), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.ϕ, ngrad[1], rtol=1e-6)
        end

        begin
            f(nθ) = dot(Δ.v, apply(CircuitGate((1,), RotationGate(nθ)), ρ).v)

            nθ = randn(3)
            ngrad = ngradient(f, nθ)
            dg = Qaintessent.backward_density(RotationGate(nθ), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)

            # special case: zero vector
            nθ = zeros(3)
            ngrad = ngradient(f, nθ)
            dg = Qaintessent.backward_density(RotationGate(nθ), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)
        end
    end

    @testset "entanglement gates" begin

        # fictitious density matrix
        ρ = DensityMatrix(randn(Float64, 16), 2)

        # fictitious gradients of cost function with respect to output density matrix
        Δ = DensityMatrix(randn(Float64, 16), 2)

        for g in [EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
            f(θ) = dot(Δ.v, apply(CircuitGate((1, 2), g(θ[])), ρ).v)
            θ = 2π*rand()
            ngrad = ngradient(f, [θ])
            dg = Qaintessent.backward_density(g(θ), reshape(kron(ρ.v, Δ.v), 16, 16))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)
        end
    end
end
