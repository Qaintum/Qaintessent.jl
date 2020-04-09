
include("../src/Qaintessent.jl")
using .Qaintessent

N = 5
cgc = CircuitGateChain{N}([
    single_qubit_circuit_gate(3, HadamardGate(), N),
    controlled_circuit_gate((1, 4), 2, RxGate(√0.2), N),
    controlled_circuit_gate((2,4), (1,5), SwapGate(), N),
    single_qubit_circuit_gate(2, PhaseShiftGate(0.2π), N),
    single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
    single_qubit_circuit_gate(1, RyGate(1.4π), N),
    two_qubit_circuit_gate(1,2, SwapGate(), N),
    controlled_circuit_gate(4, 5, TGate(), N),
    single_qubit_circuit_gate(3, SGate(), N),
])

print(cgc)
