using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent

if CUDA.functional()
    tol = 5e-2
else
    tol = 1e-5
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

@testset ExtendedTestSet "gradients for statevector" begin
    @testset ExtendedTestSet "statevector circuit gradients" begin

        # construct parametrized circuit
        N = 4
        A = rand(ComplexQ, 2 ,2)
        U, R = qr(A)
        U = Array(U)
        i = rand(1:N)

        cgc(θ, ϕ, χ, κ, ωn, η, ι, λ) = CircuitGate[
            circuit_gate(3, HadamardGate()),
            circuit_gate(2, RzGate(θ), (1, 4)),
            circuit_gate(2, TGate(), 4),
            circuit_gate(2, 3, SwapGate()),
            circuit_gate(2, SdagGate(), 1),
            circuit_gate(3, PhaseShiftGate(ϕ)),
            circuit_gate(1, RyGate(χ)),
            circuit_gate(4, HadamardGate(), 2),
            circuit_gate(3, X, 1),
            circuit_gate(2, Z, 4),
            circuit_gate(4, Y, 1),
            circuit_gate(3, 1, EntanglementYYGate(κ)),
            circuit_gate(3, RotationGate(ωn)),
            circuit_gate(3, RxGate(η), (1, 4)),
            circuit_gate(2, 3, EntanglementXXGate(ι)),
            circuit_gate(4, 1, EntanglementZZGate(λ)),
            circuit_gate(4, SGate()),
            circuit_gate(2, TdagGate(), 3),
            circuit_gate(i, MatrixGate(U)),
            circuit_gate(2, Z),
            circuit_gate(3, X),
            circuit_gate(2, Y, 4),
            circuit_gate(1, SGate()),
            circuit_gate(2, 3, SwapGate(), 4),
            circuit_gate(3, SdagGate()),
            circuit_gate(2, TGate()),
        ]

        # measurement operators
        meas(M) = MeasurementOperator.([Matrix{FloatQ}(I, 2^N, 2^N), Hermitian(M)], (Tuple(1:N),))

        # parameter values
        θ = 1.5π
        ϕ = 0.3
        χ = √2
        κ = exp(-0.4)
        η = √3
        ι = 3.33π
        λ = exp(-0.11)
        n = randn(FloatQ, 3)
        n /= norm(n)
        ωn = 0.2π * n
        M = randn(ComplexQ, 2^N, 2^N)
        M = 0.5*(M + adjoint(M))

        c = Circuit{N}(cgc(θ, ϕ, χ, κ, ωn, η, ι, λ), meas(M))
        # input quantum state
        ψ = Statevector(randn(ComplexQ, 2^N))
        ψl = deepcopy(ψ)
        # fictitious gradients of cost function with respect to circuit output
        Δ = [0.3, -1.2]

        dc = Qaintessent.gradients(c, ψ, Δ)[1]

        f(rθ, rϕ, rχ, rκ, ωn, rη, rι, rλ, M) = dot(Δ, apply!(copy!(ψl, ψ), Circuit{N}(cgc(rθ[], rϕ[], rχ[], rκ[], ωn, rη[], rι[], rλ[]), meas(M))))
        # numeric gradients
        ngrad = ngradient(f, [θ], [ϕ], [χ], [κ], ωn, [η], [ι], [λ], M)
        # symmetrize gradient with respect to measurement operator

        ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
        @test all(isapprox.(ngrad,
            (dc.moments[1][2].gate.U.θ,
            dc.moments[4][2].gate.ϕ,
            dc.moments[5][1].gate.θ,
            dc.moments[8][1].gate.θ,
            dc.moments[9][1].gate.nθ,
            dc.moments[10][1].gate.U.θ,
            dc.moments[11][1].gate.θ,
            dc.moments[11][2].gate.θ,
            sparse_matrix(dc.meas[2])), rtol=tol, atol=tol))
    end


    ##==----------------------------------------------------------------------------------------------------------------------


    @testset ExtendedTestSet "statevector circuit gradients with moments" begin

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
        meas(M) = MeasurementOperator.([Matrix{FloatQ}(I, 2^N, 2^N), Hermitian(M)], (Tuple(1:N),))

        # parameter values
        θ = 1.5π
        ϕ = 0.3
        χ = √2
        n = randn(FloatQ, 3)
        n /= norm(n)
        ωn = 0.2π * n
        M = randn(ComplexQ, 2^N, 2^N)
        M = 0.5*(M + adjoint(M))
        
        c = Circuit{N}(cgc(θ, ϕ, χ, ωn), meas(M))

        # input quantum state
        ψ = Statevector(randn(ComplexQ, 2^N))
        ψl = deepcopy(ψ)
        # fictitious gradients of cost function with respect to circuit output
        Δ = [0.3, -1.2]

        dc = Qaintessent.gradients(c, ψ, Δ)[1]

        f(rθ, rϕ, rχ, ωn, M) = dot(Δ, apply!(copy!(ψl, ψ), Circuit{N}(cgc(rθ[], rϕ[], rχ[], ωn), meas(M))))
        # numeric gradients
        ngrad = ngradient(f, [θ], [ϕ], [χ], ωn, M)
        # symmetrize gradient with respect to measurement operator
        ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
        
        @test all(isapprox.(ngrad,
            (dc.moments[1][2].gate.U.θ,
                dc.moments[3][1].gate.ϕ,
                dc.moments[4][2].gate.θ,
                dc.moments[4][1].gate.nθ,
                sparse_matrix(dc.meas[2])), rtol=tol, atol=tol))
    end
end