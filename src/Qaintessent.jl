module Qaintessent

using LinearAlgebra
using SparseArrays


include("gates.jl")
export
    AbstractGate,
    X,
    Y,
    Z,
    IGate,
    IIGate,
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

include("commute.jl")
    export
        comm

include("apply.jl")
export
    apply

include("models.jl")
export
    qft_circuit


include("gradients.jl")
include("view.jl")

include("graphs.jl")
export
    Dag,
    lshift,
    rshift,
    insert,
    insert!,
    remove,
    remove!,
    append!,
    opt_hadamard,
    opt_adjoint,
    get_controls


include("contraction_order.jl")
export
    simple_dag,
    line_graph,
    interaction_graph,
    is_tree_decomposition,
    tree_decomposition,
    contraction_order

include("python_contraction_order.jl")

end
