module Qaintessent

using LinearAlgebra
using SparseArrays

include("gates.jl")
include("circuit.jl")
include("models.jl")

export
    X,
    Y,
    Z,
    XGate,
    YGate,
    ZGate,
    HadamardGate,
    SGate,
    TGate,
    SdagGate,
    TdagGate,
    RxGate,
    RyGate,
    RzGate,
    RotationGate,
    PhaseShiftGate,
    SwapGate,
    ControlledGate,
    controlled_not,
    CircuitGate,
    single_qubit_circuit_gate,
    two_qubit_circuit_gate,
    controlled_circuit_gate,
    apply,
    CircuitGateChain,
    MeasurementOps,
    Circuit,
    qft_circuit

end
