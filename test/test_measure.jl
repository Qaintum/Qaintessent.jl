using Qaintessent: measure
using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using StatsBase


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "measurement" begin
    @testset "Measure single instance" begin
        state_a = ComplexF64[1;0;0;0]
        @test all(measure(state_a) .≈ [0;0])
    end
    
    @testset "Measure multiple instances" begin
        state_a = ComplexF64[1;0;0;0]
        circ = [circuit_gate(1, HadamardGate()), circuit_gate(2, HadamardGate())]
        state_b = apply(state_a, circ)
        ave = measure(state_b, 20000)
        @test all(isapprox(ave, [0.5;0.5], atol=1e-2))
    end
end