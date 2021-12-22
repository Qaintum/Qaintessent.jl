using Qaintessent: measure
using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using StatsBase


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "measurement" begin
    @testset "Measure single instance" begin
        state_a = ComplexQ[1;0;0;0]
        state_b = ComplexQ[1;0;0]
        state_c = ComplexQ[0.5;0;0;0]
        @test measure(state_a) == Dict(0=>1)
        @test_throws ErrorException("Shots has to be a natural integer above 0.")  measure(state_a, 0)
        @test_throws ErrorException("Statevector is not a viable quantum state (length $(length(state_b))")  measure(state_b)
        @test_throws ErrorException("Statevector is not a viable statevector (norm is $(norm(state_c))")  measure(state_c)
    end
    
    @testset "Measure multiple instances" begin
        state_a = ComplexQ[1;0;0;0]
        circ = [circuit_gate(1, HadamardGate())]
        state_b = apply(state_a, circ)
        samples = measure(state_b, 20000)
        # Measurement result should follow binomial distribution, this should return false ~1 out of a trillion
        @test (samples[0]<10500 && samples[0]>9500)
    end
end
