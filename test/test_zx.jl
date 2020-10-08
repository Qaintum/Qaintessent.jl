using Qaintessent
using BenchmarkTools

@testset ExtendedTestSet "ZX Implementation" begin
    @testset "Single Gate Construction Pauli Gates" begin
        N = 3
        cg = single_qubit_circuit_gate(1, Y, N)
        zx = ZXCircuit(cg)

        @test zx[1][].out[1][] === zx[2][]
        @test zx[1][] === zx[2][].in[1][]
    end

    @testset "Single Gate Construction Rotation Gates" begin
        N = 3
        cg = single_qubit_circuit_gate(1, RyGate(0.2), N)
        println(@btime zx = ZXCircuit($cg))
        zx = ZXCircuit(cg)

        @test zx[1][].out[1][] === zx[2][]
        @test zx[1][] === zx[2][].in[1][]
        @test zx[2][].out[1][] === zx[3][]
        @test zx[2][] === zx[3][].in[1][]
    end
end
