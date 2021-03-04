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
    apply

include("gradients.jl")
include("gradients_density.jl")

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
export unitary2circuit

end
