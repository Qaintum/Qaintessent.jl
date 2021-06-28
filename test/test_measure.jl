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
        @test measure(state_a) == Dict(0=>1)
    end
    
    @testset "Measure multiple instances" begin
        state_a = ComplexF64[1;0;0;0]
        circ = [circuit_gate(1, HadamardGate())]
        state_b = apply(state_a, circ)
        samples = measure(state_b, 20000)
        # Measurement result should follow binomial distribution, this should return false ~1 out of a trillion
        @test (samples[0]<10500 && samples[0]>9500)
    end
end
