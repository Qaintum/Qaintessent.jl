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

@testset ExtendedTestSet "gradients for density matrix" begin
    @testset ExtendedTestSet "density matrix gate gradients" begin

        @testset "single qubit gates" begin

            # fictitious density matrix
            ρ = DensityMatrix(randn(FloatQ, 4), 1)

            # fictitious gradients of cost function with respect to output density matrix
            Δ = DensityMatrix(randn(FloatQ, 4), 1)
            
            for g in [RxGate, RyGate, RzGate]
                θval = 2π*rand()

                f = (θ)->dot(Δ.v, apply(ρ, CircuitGate((1,), g(θ[]))).v)
                ngrad = ngradient(f, [θval])
                dg = Qaintessent.backward_density(g(θval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.θ, ngrad[1], rtol=tol, atol=tol)

                fa = (θ)->dot(Δ.v, Qaintessent.apply_mixed_add(ρ, CircuitGate((1,), g(θ[]))).v)
                ngrad = ngradient(fa, [θval])
                dg = Qaintessent.backward_density_mixed_add(g(θval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.θ, ngrad[1], rtol=tol, atol=tol)

                fs = (θ)->dot(Δ.v, Qaintessent.apply_mixed_sub(ρ, CircuitGate((1,), g(θ[]))).v)
                ngrad = ngradient(fs, [θval])
                dg = Qaintessent.backward_density_mixed_sub(g(θval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.θ, ngrad[1], rtol=tol, atol=tol)
            end

            for g in [XGate, YGate, ZGate, HadamardGate, SGate, SdagGate, TGate, TdagGate]
                dg = Qaintessent.backward_density(g(), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg, g())
            end

            begin
                nθval = randn(3)

                f = (nθ)->dot(Δ.v, apply(ρ, CircuitGate((1,), RotationGate(nθ))).v)
                ngrad = ngradient(f, nθval)
                dg = Qaintessent.backward_density(RotationGate(nθval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.nθ, ngrad[1], rtol=tol, atol=tol)

                fa = (nθ)->dot(Δ.v, Qaintessent.apply_mixed_add(ρ, CircuitGate((1,), RotationGate(nθ))).v)
                ngrad = ngradient(fa, nθval)
                dg = Qaintessent.backward_density_mixed_add(RotationGate(nθval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.nθ, ngrad[1], rtol=tol, atol=tol)

                fs = (nθ)->dot(Δ.v, Qaintessent.apply_mixed_sub(ρ, CircuitGate((1,), RotationGate(nθ))).v)
                ngrad = ngradient(fs, nθval)
                dg = Qaintessent.backward_density_mixed_sub(RotationGate(nθval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.nθ, ngrad[1], rtol=tol, atol=tol)

                # special case: zero vector

                ngrad = ngradient(f, zeros(3))
                dg = Qaintessent.backward_density(RotationGate(zeros(3)), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.nθ, ngrad[1], rtol=tol, atol=tol)

                ngrad = ngradient(fa, zeros(3))
                dg = Qaintessent.backward_density_mixed_add(RotationGate(zeros(3)), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.nθ, ngrad[1], rtol=tol, atol=tol)

                ngrad = ngradient(fs, zeros(3))
                dg = Qaintessent.backward_density_mixed_sub(RotationGate(zeros(3)), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.nθ, ngrad[1], rtol=tol, atol=tol)
            end

            begin
                ϕval = 2π*rand()

                f = (ϕ)->dot(Δ.v, apply(ρ, CircuitGate((1,), PhaseShiftGate(ϕ[]))).v)
                ngrad = ngradient(f, [ϕval])
                dg = Qaintessent.backward_density(PhaseShiftGate(ϕval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.ϕ, ngrad[1], rtol=tol, atol=tol)

                fa = (ϕ)->dot(Δ.v, Qaintessent.apply_mixed_add(ρ, CircuitGate((1,), PhaseShiftGate(ϕ[]))).v)
                ngrad = ngradient(fa, [ϕval])
                dg = Qaintessent.backward_density_mixed_add(PhaseShiftGate(ϕval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.ϕ, ngrad[1], rtol=tol, atol=tol)

                fs = (ϕ)->dot(Δ.v, Qaintessent.apply_mixed_sub(ρ, CircuitGate((1,), PhaseShiftGate(ϕ[]))).v)
                ngrad = ngradient(fs, [ϕval])
                dg = Qaintessent.backward_density_mixed_sub(PhaseShiftGate(ϕval), reshape(kron(ρ.v, Δ.v), 4, 4))
                @test isapprox(dg.ϕ, ngrad[1], rtol=tol, atol=tol)
            end
        end

        @testset "entanglement gates" begin

            # fictitious density matrix
            ρ = DensityMatrix(randn(FloatQ, 16), 2)

            # fictitious gradients of cost function with respect to output density matrix
            Δ = DensityMatrix(randn(FloatQ, 16), 2)

            for g in [EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
                θval = 2π*rand()

                f = (θ)->dot(Δ.v, apply(ρ, CircuitGate((1, 2), g(θ[]))).v)
                ngrad = ngradient(f, [θval])
                dg = Qaintessent.backward_density(g(θval), reshape(kron(ρ.v, Δ.v), 16, 16))
                @test isapprox(dg.θ, ngrad[1], rtol=tol, atol=tol)

                fa = (θ)->dot(Δ.v, Qaintessent.apply_mixed_add(ρ, CircuitGate((1, 2), g(θ[]))).v)
                ngrad = ngradient(fa, [θval])
                dg = Qaintessent.backward_density_mixed_add(g(θval), reshape(kron(ρ.v, Δ.v), 16, 16))
                @test isapprox(dg.θ, ngrad[1], rtol=tol, atol=tol)

                fs = (θ)->dot(Δ.v, Qaintessent.apply_mixed_sub(ρ, CircuitGate((1, 2), g(θ[]))).v)
                ngrad = ngradient(fs, [θval])
                dg = Qaintessent.backward_density_mixed_sub(g(θval), reshape(kron(ρ.v, Δ.v), 16, 16))
                @test isapprox(dg.θ, ngrad[1], rtol=tol, atol=tol)
            end
        end

        @testset "controlled gates 1" begin

            # fictitious density matrix
            ρ = DensityMatrix(randn(FloatQ, 64), 3)

            # fictitious gradients of cost function with respect to output density matrix
            Δ = DensityMatrix(randn(FloatQ, 64), 3)

            θval = 2π*rand()
            
            for g in [RxGate, RyGate, RzGate]
                
                f1(θ) = dot(Δ.v, apply(ρ, CircuitGate((1, 2, 3), ControlledGate{g}(g(θ[]), 2))).v)
                ngrad = ngradient(f1, [θval])
                dg = Qaintessent.backward_density(ControlledGate{g}(g(θval), 2), reshape(kron(ρ.v, Δ.v), 64, 64))
                #@code_warntype(Qaintessent.backward_density(ControlledGate{RyGate}(RyGate(θval), 2), reshape(kron(ρ.v, Δ.v), 64, 64)))
                @test isapprox(dg.U.θ, ngrad[1], rtol=tol, atol=tol)
                if !isapprox(dg.U.θ, ngrad[1], rtol=tol, atol=tol)
                    println(θval)
                    println(norm(dg.U.θ - ngrad[1]) / norm(dg.U.θ)) 
                end
            end

            for g in [EntanglementXXGate, EntanglementYYGate, EntanglementZZGate]
                f2(θ) = dot(Δ.v, apply(ρ, CircuitGate((1, 2, 3), ControlledGate{g}(g(θ[]), 1))).v)
                ngrad = ngradient(f2, [θval])
                dg = Qaintessent.backward_density(ControlledGate{g}(g(θval), 1), reshape(kron(ρ.v, Δ.v), 64, 64))
                @test isapprox(dg.U.θ, ngrad[1], rtol=tol, atol=tol)
                if !isapprox(dg.U.θ, ngrad[1], rtol=tol, atol=tol)
                    println(θval)
                    println(norm(dg.U.θ - ngrad[1]) / norm(dg.U.θ))
                end
            end            

            for g in [XGate, YGate, ZGate, HadamardGate, SGate, SdagGate, TGate, TdagGate]
                dg = Qaintessent.backward_density(ControlledGate{g}(g(), 2), reshape(kron(ρ.v, Δ.v), 64, 64))
                @test isapprox(dg, ControlledGate{g}(g(), 2))
            end
        end

        @testset "controlled gates 2" begin

            # fictitious density matrix
            ρ = DensityMatrix(randn(FloatQ, 256), 4)

            # fictitious gradients of cost function with respect to output density matrix
            Δ = DensityMatrix(randn(FloatQ, 256), 4)

            θval = 2π*rand()

            f1(θ) = dot(Δ.v, apply(ρ, CircuitGate((1, 2, 3, 4), ControlledGate{PhaseShiftGate}(PhaseShiftGate(θ[]), 3))).v)
            ngrad = ngradient(f1, [θval])
            dg = Qaintessent.backward_density(ControlledGate{PhaseShiftGate}(PhaseShiftGate(θval), 3), reshape(kron(ρ.v, Δ.v), 256, 256))
            @test isapprox(dg.U.ϕ, ngrad[1], rtol=tol, atol=tol)

            f2(θ) = dot(Δ.v, apply(ρ, CircuitGate((1, 2, 3, 4), ControlledGate{EntanglementYYGate}(EntanglementYYGate(θ[]), 2))).v)
            ngrad = ngradient(f2, [θval])
            dg = Qaintessent.backward_density(ControlledGate{EntanglementYYGate}(EntanglementYYGate(θval), 2), reshape(kron(ρ.v, Δ.v), 256, 256))
            @test isapprox(dg.U.θ, ngrad[1], rtol=tol, atol=tol)
        end
    end


    ##==----------------------------------------------------------------------------------------------------------------------


    @testset ExtendedTestSet "density matrix circuit gradients" begin

        @testset "individual circuit gate" begin
            # fictitious density matrix
            ρ = DensityMatrix(randn(FloatQ, 64), 3)

            # fictitious gradients of cost function with respect to output density matrix
            Δ = DensityMatrix(randn(FloatQ, 64), 3)

            θval = 2π*rand()

            cg(θ) = CircuitGate((3, 1), EntanglementYYGate(θ))
            f(θ) = dot(Δ.v, apply(ρ, cg(θ[])).v)
            ngrad = ngradient(f, [θval])
            dcg = Qaintessent.backward_density(cg(θval), ρ, Δ)
            @test isapprox(dcg.gate.θ, ngrad[1], rtol=tol, atol=tol)
        end

        @testset "moments" begin
            # fictitious density matrix
            ρ = DensityMatrix(randn(FloatQ, 256), 4)

            # fictitious gradients of cost function with respect to output density matrix
            Δ = DensityMatrix(0.1*randn(FloatQ, 256), 4)

            moments(θ, κ, ϕ) = [Moment([circuit_gate(2, RyGate(θ)), circuit_gate(3, 1, EntanglementZZGate(κ))]),
                                Moment(circuit_gate(3, PhaseShiftGate(ϕ), (1, 4)))]

            θval = 2π*rand()
            κval = 2π*rand()
            ϕval = 2π*rand()

            f(rθ, rκ, rϕ) = dot(Δ.v, apply(ρ, moments(rθ[], rκ[], rϕ[])).v)
            ngrad = ngradient(f, [θval], [κval], [ϕval])
            ρs = apply(ρ, moments(θval, κval, ϕval))
            dm, Δ0 = Qaintessent.backward_density(moments(θval, κval, ϕval), ρs, Δ)

            @test isapprox(dm[1].gates[1].gate.θ, ngrad[1], rtol=tol, atol=tol)
            @test isapprox(dm[1].gates[2].gate.θ, ngrad[2], rtol=tol, atol=tol)
            @test isapprox(dm[2].gates[1].gate.U.ϕ, ngrad[3], rtol=tol, atol=tol)
        end

        @testset "full circuit" begin

            N = 4
            A = randn(ComplexQ, 2 ,2)
            U, R = qr(A)
            U = Array(U)
            i = rand(1:N)
        
            cgc(θ, ϕ, χ, κ, ωn) = CircuitGate[
                circuit_gate(3, HadamardGate()),
                circuit_gate(2, RzGate(θ), (1, 4)),
                circuit_gate(2, TGate()),
                circuit_gate(2, 3, SwapGate()),
                circuit_gate(2, SdagGate(), 1),
                circuit_gate(3, PhaseShiftGate(ϕ)),
                circuit_gate(1, RyGate(χ)),
                circuit_gate(4, HadamardGate(), 2),
                circuit_gate(3, X, 1),
                circuit_gate(2, Z, 4),
                circuit_gate(4, Y),
                circuit_gate(3, 1, EntanglementYYGate(κ)),
                circuit_gate(3, RotationGate(ωn)),
                circuit_gate(4, SGate()),
                circuit_gate(2, TdagGate()),
                circuit_gate(i, MatrixGate(U)),
            ]
            # measurement operators
            meas(M) = MeasurementOperator.([Matrix{FloatQ}(I, 2^N, 2^N), Hermitian(M)], (Tuple(1:N),))

            # parameter values
            θval = 1.5π
            ϕval  = 0.3
            χval  = √2
            κval  = exp(-0.4)
            ωnval = 0.1π * randn(FloatQ, 3)
            Mval = randn(ComplexQ, 2^N, 2^N)
            Mval = 0.5*(Mval + adjoint(Mval))

            c = Circuit{N}(cgc(θval, ϕval, χval, κval, ωnval), meas(Mval))

            # input density matrix
            ρ = DensityMatrix(randn(FloatQ, 4^N), N)
            # fictitious gradients of cost function with respect to circuit output
            Δ = [0.3, -0.2]

            f(rθ, rϕ, rχ, rκ, ωn, M) = dot(Δ, apply(ρ, Circuit{N}(cgc(rθ[], rϕ[], rχ[], rκ[], ωn), meas(M))))
            # numeric gradients
            ngrad = ngradient(f, [θval], [ϕval], [χval], [κval], ωnval, Mval)
            # symmetrize gradient with respect to measurement operator
            ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
            dc = Qaintessent.gradients(c, ρ, Δ)[1]

            @test all(isapprox.(ngrad,
            (dc.moments[1][2].gate.U.θ,
            dc.moments[4][2].gate.ϕ,
            dc.moments[5][1].gate.θ,
            dc.moments[7][2].gate.θ,
            dc.moments[8][1].gate.nθ,
            sparse_matrix(dc.meas[2])), rtol=tol, atol=tol))
        end

        @testset "full circuit for apply!" begin

            N = 4
            A = randn(ComplexQ, 2 ,2)
            U, R = qr(A)
            U = Array(U)
            i = rand(1:N)
        
            cgc(θ, ϕ, χ, κ, ωn) = CircuitGate[
                circuit_gate(3, HadamardGate()),
                circuit_gate(2, RzGate(θ), (1, 4)),
                circuit_gate(2, TGate()),
                circuit_gate(2, 3, SwapGate()),
                circuit_gate(2, SdagGate(), 1),
                circuit_gate(3, PhaseShiftGate(ϕ)),
                circuit_gate(1, RyGate(χ)),
                circuit_gate(4, HadamardGate(), 2),
                circuit_gate(3, X, 1),
                circuit_gate(2, Z, 4),
                circuit_gate(4, Y),
                circuit_gate(3, 1, EntanglementYYGate(κ)),
                circuit_gate(3, RotationGate(ωn)),
                circuit_gate(4, SGate()),
                circuit_gate(2, TdagGate()),
                circuit_gate(i, MatrixGate(U)),
            ]
            # measurement operators
            meas(M) = MeasurementOperator.([Matrix{FloatQ}(I, 2^N, 2^N), Hermitian(M)], (Tuple(1:N),))
        
            # parameter values
            θval = 1.5π
            ϕval  = 0.3
            χval  = √2
            κval  = exp(-0.4)
            ωnval = 0.1π * randn(FloatQ, 3)
            Mval = randn(ComplexQ, 2^N, 2^N)
            Mval = 0.5*(Mval + adjoint(Mval))
        
            c = Circuit{N}(cgc(θval, ϕval, χval, κval, ωnval), meas(Mval))
        
            # input density matrix
            ρ = DensityMatrix(randn(FloatQ, 4^N), N)
            # fictitious gradients of cost function with respect to circuit output
            Δ = [0.3, -0.2]
        
            f(rθ, rϕ, rχ, rκ, ωn, M) = dot(Δ, apply(ρ, Circuit{N}(cgc(rθ[], rϕ[], rχ[], rκ[], ωn), meas(M))))
            # numeric gradients
            ngrad = ngradient(f, [θval], [ϕval], [χval], [κval], ωnval, Mval)
            # symmetrize gradient with respect to measurement operator
            ngrad[end][:] = 0.5*(ngrad[end] + adjoint(ngrad[end]))
            dc = Qaintessent.gradients!(c, apply!(ρ, c.moments), Δ)[1]
        
            @test all(isapprox.(ngrad,
            (dc.moments[1][2].gate.U.θ,
            dc.moments[4][2].gate.ϕ,
            dc.moments[5][1].gate.θ,
            dc.moments[7][2].gate.θ,
            dc.moments[8][1].gate.nθ,
            sparse_matrix(dc.meas[2])), rtol=tol, atol=tol))
        end
    end
end