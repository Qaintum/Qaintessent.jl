using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "circuits" begin
    @testset ExtendedTestSet "circuit construction" begin
        N = 3
        cgs = [
            circuit_gate(1, X),
            circuit_gate(2, HadamardGate()),
            circuit_gate(3, Z),
            circuit_gate(1, Y),
            ]
        ref_moments = Moment[
            Moment([circuit_gate(1, X),
                    circuit_gate(2, HadamardGate()),
                    circuit_gate(3, Z)]),
            Moment([circuit_gate(1, Y)])
        ]
        iwires = [(1,), (2,), (3,)]

        @testset "basic circuit with circuit gate" begin
            c = Circuit{N}(cgs[1])
            
            @test c[1] ≈ Moment(cgs[1])

            ψ = randn(ComplexQ, 2^N)
            s = Statevector(ψ)
            @test distribution(ψ, c) ≈ apply(ψ, c.moments)

            @test_throws ErrorException("Circuit does not contain any measurement operators") apply(ψ, c)

            meas = mop.([X, X, X], iwires)
            
            c = Circuit{N}(cgs[1], meas)
            ψs = apply(ψ, c.moments)

            @test c[1] ≈ Moment(cgs[1])
            @test [dot(ψs, apply(ψs, m)) for m in meas] ≈ apply(ψ, c)
            @test [dot(ψs, apply(ψs, m)) for m in meas] ≈ apply!(s, c)
        end
        
        @testset "basic circuit construction" begin
            c = Circuit{N}(cgs)
            
            for i in 1:length(c)
                @test c[i] ≈ ref_moments[i]
            end

            ψ = randn(ComplexQ, 2^N)
            s = Statevector(ψ)
            @test distribution(ψ, c) ≈ apply(ψ, c.moments)

            @test_throws ErrorException("Circuit does not contain any measurement operators") apply(ψ, c)
            @test_throws ErrorException("Circuit does not contain any measurement operators") apply!(s, c)
        end

        @testset "circuit construction without gates" begin
            meas = mop.([X, X, X], iwires)
            
            c = Circuit{N}()

            ψ = randn(ComplexQ, 2^N)
            s = Statevector(ψ)

            @test_throws ErrorException("Circuit does not contain any gates") apply(ψ, c)
            @test_throws ErrorException("Circuit does not contain any gates") apply!(s, c)

            c = Circuit{N}(meas)

            @test_throws ErrorException("Circuit does not contain any gates") apply(ψ, c)
            @test_throws ErrorException("Circuit does not contain any gates") apply!(s, c)
        end

        @testset "circuit construction with measurements" begin
            
            meas = mop.([X, X, X], iwires)
            c = Circuit{N}(cgs, meas)
            ψ = randn(ComplexQ, 2^N)
            s = Statevector(ψ)
            
            @test distribution(ψ, c) ≈ apply(ψ, c.moments)

            ψs = apply(ψ, cgs)
            @test [dot(ψs, apply(ψs, m)) for m in meas] ≈ apply(ψ, c)
            @test [dot(ψs, apply(ψs, m)) for m in meas] ≈ apply!(s, c)
        end
    end

    @testset ExtendedTestSet "circuit helper functions" begin
        N = 3
        cgs = [
            circuit_gate(1, X),
            circuit_gate(2, HadamardGate()),
            circuit_gate(3, Z),
            circuit_gate(1, Y),
            ]
        ref_moments = Moment[
            Moment([circuit_gate(1, X),
                    circuit_gate(2, HadamardGate()),
                    circuit_gate(3, Z)]),
            Moment([circuit_gate(1, Y)])
        ]
        @testset "circuit add measurement operator" begin
            c = Circuit{N}(cgs)

            meas = mop(X, 1)
            add_measurement!(c, meas)

            ψ = randn(ComplexQ, 2^N)
            s = Statevector(ψ)
            ψs = apply(ψ, cgs)

            @test [dot(ψs, apply(ψs, meas))] ≈ apply(ψ, c)
            @test [dot(ψs, apply(ψs, meas))] ≈ apply!(s, c)
        end

        @testset "circuit add measurement operators" begin
            c = Circuit{N}(cgs)

            meas = mop.([X, X, X], (1,2,3))
            add_measurement!(c, meas)

            ψ = randn(ComplexQ, 2^N)
            s = Statevector(ψ)
            ψs = apply(ψ, cgs)

            @test [dot(ψs, apply(ψs,m)) for m in meas] ≈ apply(ψ, c)
            @test [dot(ψs, apply(ψs,m)) for m in meas] ≈ apply!(s, c)
        end

        @testset "circuit add measurement operator exceptions" begin
            meas = mop.([X, X, X], (1,2,3))
            c = Circuit{N}(cgs, meas)

            meas2 = mop(Y, 1)
            @test_throws ErrorException("Measurement operators do not commute") add_measurement!(c, meas2)
            @test_throws ErrorException("Measurement operators do not commute") add_measurement!(c, [meas2])

            meas3 = mop(X, 5)
            @test_throws ErrorException("Measurement operators affecting 5 wires provided for circuit of size 3") add_measurement!(c, meas3)
            @test_throws ErrorException("Measurement operators affecting 5 wires provided for circuit of size 3") add_measurement!(c, [meas3])
        end

        @testset "circuit array operations" begin
            m1 = Moment([circuit_gate(3, ZGate(), 1), circuit_gate(2, SGate())])
            m2 = Moment([circuit_gate(2, ZGate()), circuit_gate(1, TdagGate())])
            cgc = [circuit_gate(3, ZGate(), 1), circuit_gate(2, SGate()), circuit_gate(2, ZGate()), circuit_gate(1, TdagGate())]
            c = Circuit{3}(cgc)

            @test all(c[:] .≈ c.moments)
            @test all(c[1:2] .≈ [m1, m2])
            @test c[1:2] == c.moments
            @test isempty(c) == false
            @test m2 ≈ pop!(c)
            @test isempty(c) == false
            @test m1 ≈ pop!(c)
            @test isempty(c) == true
            @test_throws ArgumentError("array must be non-empty") pop!(c)
        end

        @testset "circuit reverse" begin
            meas = mop.([X, X, X], (1,2,3))
            c = Circuit{N}(cgs, meas)
            reverse_c = reverse(c)
            reverse_ref_moments = reverse(reverse.(ref_moments))
            for i in 1:length(c)
                @test reverse_c[i] ≈ reverse_ref_moments[i]
                @test c[i] ≈ ref_moments[i]
            end

            reverse_c = reverse!(c)
            for i in 1:length(c)
                @test reverse_c[i] ≈ reverse_ref_moments[i]
                @test c[i] ≈ reverse_ref_moments[i]
            end
        end
    end
end