using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "circuit construction" begin
    N = 1
    cgs = [
        circuit_gate(1, X),
        circuit_gate(1, HadamardGate()),
        circuit_gate(1, Z),
        circuit_gate(1, Y),
        ]
    ref_moments = Moment.(cgs)
    meas = MeasurementOperator(Matrix{Float64}(I, 2^N, 2^N), (1,))
    c = Circuit{N}(cgs, [meas])
    for i in 1:length(cgs)
        @test c[i] ≈ ref_moments[i]
    end
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.moments, ψ)

    ψs = apply(cgs, ψ)
    ψref = apply(c, ψ)
    
    @test [dot(ψs, apply(meas, ψs))] ≈ apply(c, ψ)
end

@testset ExtendedTestSet "circuit construction with circuit gate measurements" begin
    N = 3
    cgc = [
        circuit_gate(1, X),
        circuit_gate(2, HadamardGate()),
        circuit_gate(3, Z),
        circuit_gate(1, Y),
    ]
    iwires = [(1,), (2,), (3,)]
    meas = MeasurementOperator.([X, X, X], iwires)
    c = Circuit{N}(cgc, meas)
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.moments, ψ)

    ψs = apply(cgc, ψ)
    @test [dot(ψs, apply(m, ψs)) for m in meas] ≈ apply(c, ψ)
end
