using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using Qaintessent.MaxKColSubgraphQAOA
using CUDA

if CUDA.functional()
    tol = 5e-2
else
    tol = 1e-4
end

##==----------------------------------------------------------------------------------------------------------------------


# adapted from https://github.com/FluxML/Zygote.jl/blob/master/test/gradcheck.jl
function ngradient(f, xs::AbstractArray...)
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in 1:length(x)
        δ = sqrt(eps(FloatQ))
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


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "max-k-col. QAOA gradients" begin
    @testset ExtendedTestSet "max-k-col. subgraph QAOA gate gradients" begin
        @testset "MaxKColGates gates" begin
            # fictitious gradients of cost function with respect to quantum gate
            Δ = randn(ComplexQ, 16, 16)
            κ = 2
            graph = Graph(2, [(1,2)])
            g  = MaxKColSubgraphPhaseSeparationGate
            f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[], κ, graph))))
            θ = 2π*rand()
            ngrad = ngradient(f, [θ])
            dg = Qaintessent.backward(g(θ, κ, graph), conj(Δ))
            @test isapprox(dg.γ, ngrad[1], rtol=1e-4, atol=1e-4)
        end

        @testset "RNearbyValuesMixerGate gates" begin
            κs = 2:6
            rs = [1, 1, 2, 1, 4, 5]
            for (κ, r) ∈ zip(κs, rs)
                Δ = randn(ComplexQ, Base.:^(2, κ), Base.:^(2, κ))
                g = RNearbyValuesMixerGate
                f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[], r, κ))))
                θ = 2π*rand()
                ngrad = ngradient(f, [θ])
                dg = Qaintessent.backward(g(θ, r, κ), conj(Δ))
                # sadly, with less tolerance it doesn't pass (could be improved by gradients that
                # are not just numerically estimated)
                @test isapprox(dg.β, ngrad[1], rtol=1e-4, atol=1e-4)
            end
        end

        @testset "ParityRingMixerGate gates" begin
            for κ in 2:6
                Δ = randn(ComplexQ, Base.:^(2, κ), Base.:^(2, κ))
                g = ParityRingMixerGate
                f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[], κ))))
                θ = 2π*rand()
                ngrad = ngradient(f, [θ])
                dg = Qaintessent.backward(g(θ, κ), conj(Δ))
                # sadly, with less tolerance it doesn't pass (could be improved by gradients that
                # are not just numerically estimated)
                @test isapprox(dg.β, ngrad[1], rtol=1e-4, atol=1e-4)
            end
        end

        @testset "PartitionMixerGate gates" begin
            κ = 4
            Δ = randn(ComplexQ, Base.:^(2, κ), Base.:^(2, κ))
            partition = [[(1, 2), (3, 4)], [(2, 3), (4, 1)], [(1, 2), (3, 4)]]
            g = PartitionMixerGate
            f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[], κ, partition))))
            θ = 2π*rand()
            ngrad = ngradient(f, [θ])
            dg = Qaintessent.backward(g(θ, κ, partition), conj(Δ))
            @test isapprox(dg.β, ngrad[1], rtol=1e-4, atol=1e-4)
        end
    end
end