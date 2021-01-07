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
    MatrixGate,
    matrix

include("circuitgate.jl")
export
    AbstractCircuitGate,
    CircuitGate,
    circuit_gate

include("circuit.jl")
export
    Moment,
    single_qubit_circuit_gate,
    two_qubit_circuit_gate,
    controlled_circuit_gate,
    MeasurementOperator,
    Circuit,
    distribution,
    rdm

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

include("gradients.jl")
include("view.jl")

# include("graphs.jl")
# export
#     Dag,
#     optimize

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

include("compile.jl")

end
