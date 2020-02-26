module Qaintessent

using LinearAlgebra

include("gates.jl")

export
    X,
    Y,
    Z,
    XGate,
    YGate,
    ZGate,
    SGate,
    TGate,
    SdagGate,
    TdagGate,
    SwapGate,
    ControlledGate,
    controlled_not


end
