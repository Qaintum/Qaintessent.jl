module Qaintessent

using LinearAlgebra
using SparseArrays


include("util.jl")

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
    EntanglementXXGate,
    EntanglementYYGate,
    EntanglementZZGate,
    ControlledGate,
    controlled_not,
    MatrixGate

include("register.jl")
export
    qreg,
    creg,
    reg_check,
    set_creg!,
    add_creg!

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

include("density_matrix.jl")
export
    DensityMatrix,
    pauli_group_matrix,
    density_from_statevector

include("commute.jl")
    export
        iscommuting

include("apply.jl")
export
    apply

include("models.jl")
export
    qft_circuit,
    toffoli_circuit,
    vbe_adder_circuit,
    qcla_out_adder_circuit,
    qcla_inplace_adder_circuit

include("gradients.jl")
include("view.jl")

include("graphs.jl")
export
    Dag,
    optimize!

include("openqasm/grammar.jl")
include("openqasm/gate_transformers.jl")

include("openqasm/transform_qasm.jl")
export
    qasm2cgc

include("openqasm/transform_cgc.jl")
export
    cgc2qasm

include("fock.jl")
export
    FockState,
    vacuum_fock_state,
    coherent_state,
    â,
    âDag,
    n̂,
    D̂,
    Û,
    Ŝ,
    wigner

include("compile.jl")

end
