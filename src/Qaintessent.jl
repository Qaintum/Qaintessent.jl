module Qaintessent

using LinearAlgebra
using SparseArrays
using Memoize
using CUDA

if CUDA.functional()
    const FloatQ = Float32
    const ComplexQ = ComplexF32
    to_gpu(x) = CuArray(x)
    promote(f::Float64) = convert(Float32, f)
    promote(c::ComplexF64) = convert(ComplexF32, c)
else
    const FloatQ = Float64
    const ComplexQ = ComplexF64
    to_gpu(x) = x
    promote(f::Float64) = convert(Float64, f)
    promote(c::ComplexF64) = convert(ComplexF64, c)
end

export FloatQ, ComplexQ


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
    num_wires, 
    data

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

module gpu
    include("gpu/apply_gpu.jl")
end

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
        Graph,
        max_k_col_subgraph_phase_separation_hamiltonian

    include("qaoa/qaoa_gradients.jl")
end

end
