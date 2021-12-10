module Qaintessent

using LinearAlgebra
using SparseArrays
using Memoize


include("util.jl")
export 
    rdm

include("gates.jl")
export
    AbstractGate,
    X,
    Y,
    Z,
    H,
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
    EntanglementXXGate,
    EntanglementYYGate,
    EntanglementZZGate,
    ControlledGate,
    controlled_not,
    MatrixGate,
    matrix,
    sparse_matrix,
    num_wires

include("circuitgate.jl")
export
    AbstractCircuitGate,
    CircuitGate,
    circuit_gate,
    req_wires

include("moment.jl")
export
    Moment

include("measurementoperator.jl")
export
    MeasurementOperator,
    mop

include("circuit.jl")
export
    Circuit,
    add_measurement!,
    distribution

include("statevector.jl")
export 
    Statevector

include("density_matrix.jl")
export
    DensityMatrix,
    pauli_group_matrix,
    density_from_statevector,
    density_from_matrix

include("commute.jl")
    export
        iscommuting

include("apply.jl")
export
    apply
    
include("apply_density.jl")
export
    apply,
    apply!

include("apply_statevector.jl")
export
    apply!

include("gradients.jl")
include("gradients_density.jl")
include("gradients_statevector.jl")

include("view.jl")

include("openqasm/register.jl")
export
    qreg,
    creg,
    reg_check,
    set_creg!,
    add!,
    add_control!

include("openqasm/grammar.jl")
include("openqasm/gate_transformers.jl")

include("openqasm/transform_qasm.jl")
export
    qasm2cgc

include("openqasm/transform_cgc.jl")
export
    cgc2qasm

include("unitary2circuit.jl")
export 
    unitary2circuit

include("measure.jl")
export 
    measure


module QAOAHelperDataStructs
    include("qaoa/qaoa_helper_data_structs.jl")
    export
        Graph,
        EdgeWeightedGraph,
        to_edge_weighted_graph,
        adjacency_matrix
end

module MaxKColSubgraphQAOA
    include("qaoa/mixer_gates.jl")
    export
        RNearbyValuesMixerGate,
        ParityRingMixerGate,
        PartitionMixerGate,
        r_nearby_values_hamiltonian_onehot

    include("qaoa/phase_separator_gates.jl")
    export
        MaxKColSubgraphPhaseSeparationGate,
        max_k_col_subgraph_phase_separation_hamiltonian

    include("qaoa/qaoa_gradients.jl")
end

module MaxCutWSQAOA
    include("qaoa/phase_separator_gates.jl")
    export
        MaxCutPhaseSeparationGate,
        max_cut_phase_separation_hamiltonian
    
    include("qaoa/mixer_gates.jl")
    export
        WSQAOAMixerGate,
        wsqaoa_mixer_hamiltonian,
        RxMixerGate,
        rx_mixer_hamiltonian

    include("qaoa/qaoa_gradients.jl")
end

using .QAOAHelperDataStructs
using .MaxKColSubgraphQAOA
using .MaxCutWSQAOA

end
