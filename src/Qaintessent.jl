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
    SwapGate,
    ControlledGate,
    controlled_not,
    CircuitGate

end
