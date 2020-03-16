module Qaintessent

using LinearAlgebra
using SparseArrays

include("gates.jl")
include("circuit.jl")

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
    SwapGate,
    ControlledGate,
    controlled_not,
    CircuitGate,
    single_qubit_circuit_gate,
    two_qubit_circuit_gate,
    controlled_circuit_gate,
    CircuitBlock

end
