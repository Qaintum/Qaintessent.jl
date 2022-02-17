using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using Qaintessent.QAOAHelperDataStructs
using Qaintessent.MaxCutWSQAOA

##==----------------------------------------------------------------------------------------------------------------------


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


##==----------------------------------------------------------------------------------------------------------------------

@testset ExtendedTestSet "maxcut WS-QAOA gate gradients" begin
    @testset "MaxCutPhaseSeparationGate gates" begin
        Δ = randn(ComplexF64, 16, 16)
        graph = Graph(4, [(1, 2), (2, 3), (2, 4), (3, 4)])
        g  = MaxCutPhaseSeparationGate
        f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[], graph))))
        θ = 2π*rand()
        ngrad = ngradient(f, [θ])
        dg = Qaintessent.backward(g(θ, graph), conj(Δ))
        @test isapprox(dg.γ, ngrad[1], rtol=1e-5, atol=1e-5)
    end

    @testset "WSQAOAMixerGate gates" begin
        Δ = randn(ComplexF64, 16, 16)
        c_opt = [0.0, 0.0, 1.0, 1.0]
        g = WSQAOAMixerGate
        f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[], c_opt, 0.0, false))))
        θ = 2π*rand()
        ngrad = ngradient(f, [θ])
        dg = Qaintessent.backward(g(θ, c_opt, 0.0, false), conj(Δ))
        @test isapprox(dg.β, ngrad[1], rtol=1e-5, atol=1e-5)
    end
end