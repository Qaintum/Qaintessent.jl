module Qaintessent

using LinearAlgebra
using SparseArrays


include("gates.jl")
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
    controlled_not


include("circuit.jl")
export
    CircuitGate,
    single_qubit_circuit_gate,
    two_qubit_circuit_gate,
    controlled_circuit_gate,
    apply,
    rdm,
    CircuitGateChain,
    MeasurementOps,
    Circuit


include("models.jl")
export
    qft_circuit


include("gradients.jl")
include("view.jl")

end
