# Circuit Construction and Usage

The output of typical Quantum simulation is the expectation values of the output quantum state based on certain measurement operators. In Qaintessent.jl, this is accomplished by applying a Quantum [Circuit](@ref)(s) consisting of a  [CircuitGateChain](@ref) and [Measurement Operators](@ref). See the corresponding sections for more information.

```@meta
CurrentModule = Qaintessent
```

## CircuitGate
The basic building block of any [CircuitGateChain](@ref) is CircuitGates. These are constructed from basic quantum [Gates](@ref). The constructor for CircuitGates is complicated: [CircuitGate Helper Functions](@ref) and a [CircuitGate Example](@ref) can be found below.

```@docs
CircuitGate
CircuitGate{M,N,G}(::NTuple{M, <:Integer}, ::G) where {M,N,G}
CircuitGate(::NTuple{M, <:Integer}, ::AbstractGate{M}, N::Integer) where {M}
apply(cg::CircuitGate{M,N,G}, ψ::AbstractVector) where {M,N,G}
```

### CircuitGate Helper Functions
```@docs
single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate{1}, N::Integer)
two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate{2}, N::Integer)
controlled_circuit_gate
```

### CircuitGate Example

The `Qaintessent.matrix` function can be used to convert `CircuitGate` objects to a CSC sparse matrix representation

```@example CircuitGate
    using Qaintessent
    N = 2
    cnot = controlled_circuit_gate(1, 2, XGate(), N)

    Qaintessent.matrix(cnot)
```
The `CircuitGate` can then be applied to a quantum state in state vector form.

```@example CircuitGate
    ψ = [0, 0, 1, 0]

    apply(cnot, ψ)
```

## Moments
When performing a quantum simulation, it may be required to define an intermediate state; this can be used to optimize the mapping of gates to physical qubits. In Qaintessent.jl, we allow for groupings of [CircuitGate](@ref)(s) that can be run in parallel into [Moments](@ref). A [Moment Example](@ref) is shown below.

```@docs
Moment
Moment{N}(g::AbstractCircuitGate{N}) where {N}
apply(m::Moment{N}, ψ::AbstractVector) where {N}
```

### Moment Example

```@example Moment
    using Qaintessent

    N = 3
    m = Moment{N}(
        [
            controlled_circuit_gate(1, 2, XGate(), N),
            single_qubit_circuit_gate(3, RxGate(0.2π), N),
        ]
    )
    println(m) # hide
```
The `Moment` can then be applied to a quantum state in state vector form.

```@example Moment
    ψ = randn(ComplexF64, 2^N)
    
    apply(m, ψ)
```

## CircuitGateChain
CircuitGateChains can be constructed in two ways, either by providing a Vector of `CircuitGate` objects or `Moment` objects. A [CircuitGateChain Example](@ref) can be seen below.

```@docs
CircuitGateChain
CircuitGateChain{N}(gates::AbstractVector{<:AbstractCircuitGate{N}}) where {N}
apply(c::CircuitGateChain{N}, ψ::AbstractVector) where {N}
```

### CircuitGateChain Example
In the following example, a 3-qubit Quantum Fourier Transform is applied to an arbitrary input state

```@example CGC
    using Qaintessent

    N = 3
    cgc = CircuitGateChain{N}(
        [
            single_qubit_circuit_gate(1, HadamardGate(), N),
            controlled_circuit_gate(2, 1, PhaseShiftGate(2π/2^2), N),
            controlled_circuit_gate(3, 1, PhaseShiftGate(2π/2^3), N),
            single_qubit_circuit_gate(2, HadamardGate(), N),
            controlled_circuit_gate(3, 1, PhaseShiftGate(2π/2^2), N),
            single_qubit_circuit_gate(2, HadamardGate(), N),
        ]
    )
    println(cgc) # hide
```
The `CircuitGateChain` can then be applied to a quantum state in state vector form.

```@example CGC
    ψ = randn(ComplexF64, 2^N)
    
    apply(cgc, ψ)
```

## Measurement Operators
Measurement operators can be defined in a `MeasurementOps` object. These are used in conjunction with a [CircuitGateChain](@ref) to create a [Circuit](@ref) object. Usage of `MeasurementOps` can be seen in the [Circuit Example](@ref) shown below.

```@docs
MeasurementOps
MeasurementOps{N}(mop::AbstractMatrix) where {N}
```

## Circuit
`Circuit` objects combine a `CircuitGateChain` and a `MeasurementOps` object. Applying a `Circuit` to a given input quantum state outputs the various expectation values from the measurement operators defined in the `MeasurementOps` object. A simple circuit is shown in the [Circuit Example](@ref).

```@docs
Circuit
apply(c::Circuit{N}, ψ::AbstractVector) where {N}
```

### Circuit Example

```@example Circuit
    using Qaintessent
    N = 3
    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(1, HadamardGate(), N),   
        single_qubit_circuit_gate(2, HadamardGate(), N),
        single_qubit_circuit_gate(3, HadamardGate(), N),
    ])
    I = [1 0 ; 0 1]
    x = Qaintessent.matrix(XGate())
    meas_op = kron(kron(I,x), I)
    meas = MeasurementOps{N}([meas_op])
    c = Circuit{N}(cgc, meas)
    println(c) # hide
```
Applying the `Circuit` object to a 3-qubit quantum state all in the ground 0 state.

```@example Circuit
    ψ = Float64[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    
    apply(c, ψ)
```