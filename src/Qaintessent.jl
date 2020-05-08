module Qaintessent

using LinearAlgebra
using SparseArrays


include("gates.jl")
export
    AbstractGate,
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
    ishermitian


include("circuit.jl")
export
    AbstractCircuitGate,
    CircuitGate,
    AbstractMoment,
    Moment,
    single_qubit_circuit_gate,
    two_qubit_circuit_gate,
    controlled_circuit_gate,
    rdm,
    CircuitGateChain,
    MeasurementOps,
    Circuit,
    distribution

include("apply.jl")
export
    apply

include("models.jl")
export
    qft_circuit


include("gradients.jl")
include("view.jl")

end
