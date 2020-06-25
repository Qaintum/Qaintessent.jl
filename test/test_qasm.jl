using Test
using TestSetExtensions
using Qaintessent


@testset ExtendedTestSet "test qasm" begin
    filename = "test4.qasm"
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
        single_qubit_circuit_gate(6, RxGate(0.0), N),
        single_qubit_circuit_gate(6, RyGate(0.1), N),
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
