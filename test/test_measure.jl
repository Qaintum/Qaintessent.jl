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
        ave = measure(state_b, 20000)
    end
end
