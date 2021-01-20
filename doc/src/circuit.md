# Circuit Construction and Usage

The output of typical Quantum simulation is the expectation values of the output quantum state based on certain measurement operators. In Qaintessent.jl, this is accomplished by applying a Quantum [Circuit](@ref)(s) consisting of a  series of [CircuitGate](@ref) objects and [Measurement Operators](@ref). See the corresponding sections for more information.

```@meta
CurrentModule = Qaintessent
```

## CircuitGate
The basic building block of any [Circuit](@ref) is CircuitGates. These are constructed from basic quantum [Gates](@ref). The constructor for CircuitGates is complicated: [CircuitGate Helper Functions](@ref) and a [CircuitGate Example](@ref) can be found below.

```@docs
CircuitGate
apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}
```

### CircuitGate Helper Functions
```@docs
circuit_gate
```

### CircuitGate Example

The `Qaintessent.matrix` function can be used to convert `CircuitGate` objects to a CSC sparse matrix representation

```@example CircuitGate
    using Qaintessent
    N = 2
    cnot = circuit_gate(1, XGate(), 2)

    Qaintessent.sparse_matrix(cnot)
```
The `CircuitGate` can then be applied to a quantum state in state vector form.

```@example CircuitGate
    ψ = ComplexF64[0, 0, 1, 0]

    apply(cnot, ψ)
```

## Moments
When performing a quantum simulation, it may be required to define an intermediate state; this can be used to optimize the mapping of gates to physical qubits. In Qaintessent.jl, we allow for groupings of [CircuitGate](@ref)(s) that can be run in parallel into [Moments](@ref). A [Moment Example](@ref) is shown below.

```@docs
Moment
apply(m::Moment, ψ::Vector{<:Complex})
```

### Moment Example

```@example Moment
    using Qaintessent
    N = 3
    m = Moment(
        [
            circuit_gate(1, X, 2),
            circuit_gate(3, RxGate(0.2π)),
        ]
    )
    println(m) # hide
```
The `Moment` can then be applied to a quantum state in state vector form.

```@example Moment
    ψ = randn(ComplexF64, 2^N)

    apply(m, ψ)
```

## Measurement Operators
Measurement operators can be defined in a `MeasurementOperator` object. These are used in conjunction with a series of[CircuitGate](@ref) objects to create a [Circuit](@ref) object. Usage of `MeasurementOperator` can be seen in the [Circuit Example](@ref) shown below.

```@docs
MeasurementOperator
```

## Circuit
`Circuit` objects combine a Vector of `Moment` and a `MeasurementOperator` objects. Applying a `Circuit` to a given input quantum state outputs the various expectation values from the measurement operators defined in the `MeasurementOperator` objects. A simple circuit is shown in the [Circuit Example](@ref).

```@docs
Circuit
apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}
```

### Circuit Example

```@example Circuit
    using Qaintessent
    N = 3
    cgc = CircuitGate[
        circuit_gate(1, HadamardGate()),   
        circuit_gate(2, HadamardGate()),
        circuit_gate(3, HadamardGate()),
    ]
    I = ComplexF64[1 0 ; 0 1]
    meas = [MeasurementOperator(X, (1,))]
    c = Circuit{N}(cgc, meas)
    println(c) # hide
```
Applying the `Circuit` object to a 3-qubit quantum state all in the ground 0 state.

```@example Circuit
    ψ = ComplexF64[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    apply(c, ψ)
```
