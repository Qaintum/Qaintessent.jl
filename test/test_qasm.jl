using Test
using TestSetExtensions
using Qaintessent
using LinearAlgebra

@testset ExtendedTestSet "qasm import test" begin
    filename = "import_test.qasm"
    c = import_file(filename)
    N = 10
    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(2, XGate(), N),
        single_qubit_circuit_gate(6, XGate(), N),
        single_qubit_circuit_gate(7, XGate(), N),
        single_qubit_circuit_gate(8, XGate(), N),
        single_qubit_circuit_gate(9, XGate(), N),
        controlled_circuit_gate((1),(6), XGate(), N),
        controlled_circuit_gate((1),(2), XGate(), N),
        controlled_circuit_gate((1,2),(6), XGate(), N),
        controlled_circuit_gate((5),(10), XGate(), N),
        single_qubit_circuit_gate(4, HadamardGate(), N),
        single_qubit_circuit_gate(9, HadamardGate(), N),
        single_qubit_circuit_gate(5, HadamardGate(), N),
        single_qubit_circuit_gate(6, RzGate(0.1), N),
        single_qubit_circuit_gate(6, RyGate(0.0), N),
        single_qubit_circuit_gate(6, RzGate(0.3), N),
        controlled_circuit_gate((3),(6), XGate(), N),
        controlled_circuit_gate((5,3),(6), XGate(), N)]
        ; creg=[0,0,0,0,0,0,0])
    meas = MeasurementOps{N}(AbstractMatrix[])
    c_ref = Circuit{N}(cgc_ref, meas)
    for i in 1:length(c_ref.cgc)
        @test c_ref[i] ≈ c[i]
    end
end

@testset ExtendedTestSet "qasm adder import test" begin
    filename = "import_adder_test.qasm"
    N = 10
    c = import_file(filename)
    ψ = fill(0.0+0.0im, 2^N)
    ψ[1] = 1.0+0.0im

    result = apply(c, ψ)
    answer = 0
    
    for (i,flag) in enumerate(result)
        if flag == -1
            answer += 2^i
        end
    end
    @test answer == 32
end

@testset ExtendedTestSet "qasm export test" begin
    filename = "export.qasm"
    ref_filename = "export_ref.qasm"
    N = 4
    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, XGate(), N),
        two_qubit_circuit_gate(3, 2, SwapGate(), N),
        single_qubit_circuit_gate(1, RyGate(0.3π), N),
        single_qubit_circuit_gate(2, RzGate(2.4), N),
        single_qubit_circuit_gate(2, TGate(), N),
        controlled_circuit_gate(1,2, XGate(), N),
        controlled_circuit_gate(3,1, RzGate(0.3), N)
        ])
    meas = MeasurementOps{N}(AbstractMatrix[])
    c_ref = Circuit{N}(cgc_ref, meas)
    export_file(c_ref, filename)
    f = open(filename)
    g = open(ref_filename)

    @test String(read(f)) == String(read(g))
    close(f)
    close(g)
    rm(filename)
    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, XGate(), N),
        two_qubit_circuit_gate(3, 2, SwapGate(), N),
        controlled_circuit_gate((2,3),1, RzGate(0.3), N)
        ])
    c_ref = Circuit{N}(cgc_ref, meas)
    @test_throws ErrorException("Gate CircuitGate{3,4,ControlledGate{1,3}}((2, 3, 1), ControlledGate{1,3}(RzGate([0.3])), Int64[]) has 2 control wires. Multi-control wires are currently not supported for the OpenQASM format.") export_file(c_ref, filename)
    cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, XGate(), N),
        two_qubit_circuit_gate(3, 2, SwapGate(), N),
        controlled_circuit_gate((2),1, RyGate(0.3), N)
        ])
    c_ref = Circuit{N}(cgc_ref, meas)
    @test_throws ErrorException("Controlled Gate RyGate conversion to OpenQASM is currently not supported in Qaintessent.jl") export_file(c_ref, filename)
end
