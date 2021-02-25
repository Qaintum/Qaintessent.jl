using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


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


@testset ExtendedTestSet "gate gradients" begin

    @testset "single qubit gates" begin

        # fictitious gradients of cost function with respect to quantum gate
        Δ = randn(ComplexF64, 2, 2)

        for g in [RxGate, RyGate, RzGate]
            f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[]))))
            θ = 2π*rand()
            ngrad = ngradient(f, [θ])
            dg = Qaintessent.backward(g(θ), conj(Δ))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-5)
        end

        begin
            f(ϕ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(PhaseShiftGate(ϕ[]))))
            ϕ = 2π*rand()
            ngrad = ngradient(f, [ϕ])
            dg = Qaintessent.backward(PhaseShiftGate(ϕ), conj(Δ))
            @test isapprox(dg.ϕ, ngrad[1], rtol=1e-5)
        end

        begin
            f(nθ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(RotationGate(nθ))))
            nθ = randn(3)
            ngrad = ngradient(f, nθ)
            dg = Qaintessent.backward(RotationGate(nθ), conj(Δ))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)
        end

        begin
            f(nθ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(RotationGate(nθ))))
            # special case: zero vector
            nθ = zeros(3)
            ngrad = ngradient(f, nθ)
            dg = Qaintessent.backward(RotationGate(nθ), conj(Δ))
            @test isapprox(dg.nθ, ngrad[1], rtol=1e-6)
        end
    end

    @testset "entanglement gates" begin
        # fictitious gradients of cost function with respect to quantum gate
        Δ = randn(ComplexF64, 4, 4)

        for g in [EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
            f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(g(θ[]))))
            θ = 2π*rand()
            ngrad = ngradient(f, [θ])
            dg = Qaintessent.backward(g(θ), conj(Δ))
            @test isapprox(dg.θ, ngrad[1], rtol=1e-5)
        end
    end

    @testset "controlled gates" begin
        # fictitious gradients of cost function with respect to quantum gate
        Δ = randn(ComplexF64, 8, 8)

        f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(ControlledGate{RyGate}(RyGate(θ[]), 2))))
        θ = 2π*rand()
        ngrad = ngradient(f, [θ])
        dg = Qaintessent.backward(ControlledGate{RyGate}(RyGate(θ), 2), conj(Δ))
        @test isapprox(dg.U.θ, ngrad[1], rtol=1e-5)

        f(θ) = 2*real(sum(Δ .* Qaintessent.sparse_matrix(ControlledGate{EntanglementXXGate}(EntanglementXXGate(θ[]), 1))))
        θ = 2π*rand()
        ngrad = ngradient(f, [θ])
        dg = Qaintessent.backward(ControlledGate{EntanglementXXGate}(EntanglementXXGate(θ), 2), conj(Δ))
        @test isapprox(dg.U.θ, ngrad[1], rtol=1e-5)
    end
end


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "circuit gradients" begin

    # construct parametrized circuit
    N = 4
    A = rand(ComplexF64, 2 ,2)
    i = rand(1:N)
    U, R = qr(A)
    U = Array(U)

    cgc(θ, ϕ, χ, κ, ωn) = CircuitGate[
        circuit_gate(3, HadamardGate()),
        circuit_gate(2, RzGate(θ), (1, 4)),
        circuit_gate(2, 3, SwapGate()),
        circuit_gate(3, PhaseShiftGate(ϕ)),
        circuit_gate(3, RotationGate(ωn)),
        circuit_gate(1, RyGate(χ)),
        circuit_gate(3, 1, EntanglementYYGate(κ)),
        circuit_gate(i, MatrixGate(U))
    ]
    # measurement operators
    meas(M) = MeasurementOperator.([Matrix{Float64}(I, 2^N, 2^N), Hermitian(M)], (Tuple(1:N),))

    # parameter values
    θ = 1.5π
    ϕ = 0.3
    χ = √2
    κ = exp(-0.4)
    n = randn(Float64, 3)
    n /= norm(n)
    ωn = 0.2π * n
    M = randn(ComplexF64, 2^N, 2^N)
    M = 0.5*(M + adjoint(M))

    c = Circuit{N}(cgc(θ, ϕ, χ, κ, ωn), meas(M))
    # input quantum state
    ψ = randn(ComplexF64, 2^N)
    ψl = deepcopy(ψ)
    # fictitious gradients of cost function with respect to circuit output
    Δ = [0.3, -1.2]

    dc = Qaintessent.gradients(c, ψ, Δ)[1]

    f(rθ, rϕ, rχ, rκ, ωn, M) = dot(Δ, apply(Circuit{N}(cgc(rθ[], rϕ[], rχ[], rκ[], ωn), meas(M)), ψ))
    # numeric gradients
    ngrad = ngradient(f, [θ], [ϕ], [χ], [κ], ωn, M)
    # symmetrize gradient with respect to measurement operator
    ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))

    @test all(isapprox.(ngrad,
        (dc.moments[1][2].gate.U.θ,
         dc.moments[3][1].gate.ϕ,
         dc.moments[4][2].gate.θ,
         dc.moments[5][1].gate.θ,
         dc.moments[4][1].gate.nθ,
         sparse_matrix(dc.meas[2])), rtol=1e-5, atol=1e-5))
end


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "circuit gradients with moments" begin

   # construct parametrized circuit
   N = 4
    cgc(θ, ϕ, χ, ωn) = Moment[
        Moment([
            circuit_gate(3, HadamardGate()),
            circuit_gate(2, RzGate(θ), (1, 4)),
        ]),
        Moment(circuit_gate(2, 3, SwapGate())),
        Moment(circuit_gate(3, PhaseShiftGate(ϕ))),
        Moment([
            circuit_gate(3, RotationGate(ωn)),
            circuit_gate(1, RyGate(χ)),
        ]),
    ]
   # measurement operators
   meas(M) = MeasurementOperator.([Matrix{Float64}(I, 2^N, 2^N), Hermitian(M)], (Tuple(1:N),))

   # parameter values
   θ = 1.5π
   ϕ = 0.3
   χ = √2
   n = randn(Float64, 3)
   n /= norm(n)
   ωn = 0.2π * n
   M = randn(ComplexF64, 2^N, 2^N)
   M = 0.5*(M + adjoint(M))
   
   c = Circuit{N}(cgc(θ, ϕ, χ, ωn), meas(M))

   # input quantum state
   ψ = randn(ComplexF64, 2^N)

   # fictitious gradients of cost function with respect to circuit output
   Δ = [0.3, -1.2]

   dc = Qaintessent.gradients(c, ψ, Δ)[1]

   f(rθ, rϕ, rχ, ωn, M) = dot(Δ, apply(Circuit{N}(cgc(rθ[], rϕ[], rχ[], ωn), meas(M)), ψ))
   # numeric gradients
   ngrad = ngradient(f, [θ], [ϕ], [χ], ωn, M)
   # symmetrize gradient with respect to measurement operator
   ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
   
   @test all(isapprox.(ngrad,
       (dc.moments[1][2].gate.U.θ,
        dc.moments[3][1].gate.ϕ,
        dc.moments[4][2].gate.θ,
        dc.moments[4][1].gate.nθ,
        sparse_matrix(dc.meas[2])), rtol=1e-5, atol=1e-5))
end
