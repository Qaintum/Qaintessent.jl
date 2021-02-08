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
            θval = 2π*rand()

            f(θ) = dot(Δ.v, apply(CircuitGate((1,), g(θ[])), ρ).v)
            ngrad = ngradient(f, [θval])
            dg = Qaintessent.backward_density(g(θval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)

            fa(θ) = dot(Δ.v, Qaintessent.apply_mixed_add(CircuitGate((1,), g(θ[])), ρ).v)
            ngrad = ngradient(fa, [θval])
            dg = Qaintessent.backward_density_mixed_add(g(θval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)

            fs(θ) = dot(Δ.v, Qaintessent.apply_mixed_sub(CircuitGate((1,), g(θ[])), ρ).v)
            ngrad = ngradient(fs, [θval])
            dg = Qaintessent.backward_density_mixed_sub(g(θval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)
        end

        begin
            nθval = randn(3)

            f(nθ) = dot(Δ.v, apply(CircuitGate((1,), RotationGate(nθ)), ρ).v)
            ngrad = ngradient(f, nθval)
            dg = Qaintessent.backward_density(RotationGate(nθval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)

            fa(nθ) = dot(Δ.v, Qaintessent.apply_mixed_add(CircuitGate((1,), RotationGate(nθ)), ρ).v)
            ngrad = ngradient(fa, nθval)
            dg = Qaintessent.backward_density_mixed_add(RotationGate(nθval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)

            fs(nθ) = dot(Δ.v, Qaintessent.apply_mixed_sub(CircuitGate((1,), RotationGate(nθ)), ρ).v)
            ngrad = ngradient(fs, nθval)
            dg = Qaintessent.backward_density_mixed_sub(RotationGate(nθval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)

            # special case: zero vector

            ngrad = ngradient(f, zeros(3))
            dg = Qaintessent.backward_density(RotationGate(zeros(3)), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)

            ngrad = ngradient(fa, zeros(3))
            dg = Qaintessent.backward_density_mixed_add(RotationGate(zeros(3)), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)

            ngrad = ngradient(fs, zeros(3))
            dg = Qaintessent.backward_density_mixed_sub(RotationGate(zeros(3)), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)
        end

        begin
            ϕval = 2π*rand()

            f(ϕ) = dot(Δ.v, apply(CircuitGate((1,), PhaseShiftGate(ϕ[])), ρ).v)
            ngrad = ngradient(f, [ϕval])
            dg = Qaintessent.backward_density(PhaseShiftGate(ϕval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.ϕ, ngrad[1], rtol=1e-6)

            fa(ϕ) = dot(Δ.v, Qaintessent.apply_mixed_add(CircuitGate((1,), PhaseShiftGate(ϕ[])), ρ).v)
            ngrad = ngradient(fa, [ϕval])
            dg = Qaintessent.backward_density_mixed_add(PhaseShiftGate(ϕval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.ϕ, ngrad[1], rtol=1e-6)

            fs(ϕ) = dot(Δ.v, Qaintessent.apply_mixed_sub(CircuitGate((1,), PhaseShiftGate(ϕ[])), ρ).v)
            ngrad = ngradient(fs, [ϕval])
            dg = Qaintessent.backward_density_mixed_sub(PhaseShiftGate(ϕval), reshape(kron(ρ.v, Δ.v), 4, 4))
            @test isapprox(dg.ϕ, ngrad[1], rtol=1e-6)
        end
    end

    @testset "entanglement gates" begin

        # fictitious density matrix
        ρ = DensityMatrix(randn(Float64, 16), 2)

        # fictitious gradients of cost function with respect to output density matrix
        Δ = DensityMatrix(randn(Float64, 16), 2)

        for g in [EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
            θval = 2π*rand()

            f(θ) = dot(Δ.v, apply(CircuitGate((1, 2), g(θ[])), ρ).v)
            ngrad = ngradient(f, [θval])
            dg = Qaintessent.backward_density(g(θval), reshape(kron(ρ.v, Δ.v), 16, 16))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)

            fa(θ) = dot(Δ.v, Qaintessent.apply_mixed_add(CircuitGate((1, 2), g(θ[])), ρ).v)
            ngrad = ngradient(fa, [θval])
            dg = Qaintessent.backward_density_mixed_add(g(θval), reshape(kron(ρ.v, Δ.v), 16, 16))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)

            fs(θ) = dot(Δ.v, Qaintessent.apply_mixed_sub(CircuitGate((1, 2), g(θ[])), ρ).v)
            ngrad = ngradient(fs, [θval])
            dg = Qaintessent.backward_density_mixed_sub(g(θval), reshape(kron(ρ.v, Δ.v), 16, 16))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-6)
        end
    end

end
