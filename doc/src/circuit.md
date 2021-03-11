# Circuit Construction and Usage

The output of typical Quantum simulation is the expectation values of the output quantum state based on certain measurement operators. In Qaintessent.jl, this is accomplished by applying a Quantum [Circuit](@ref) consisting of a  series of [`CircuitGate`](@ref) objects and [Measurement Operators](@ref).

```@meta
CurrentModule = Qaintessent
DocTestSetup = quote
    using Qaintessent
end
```

## Moments
When performing a quantum simulation, it may be required to define an intermediate state; this can be used to optimize the mapping of gates to physical qubits. In Qaintessent.jl, we allow for groupings of [CircuitGate](@ref)(s) that can be run in parallel into [Moments](@ref). A [Moment Example](@ref) is shown below.

```@docs
Moment
apply(ψ::Vector{<:Complex}, m::Moment)
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

    apply(ψ, m)
```

## Measurement Operators
Measurement operators can be defined in a `MeasurementOperator` object. These are used in conjunction with a series of[CircuitGate](@ref) objects to create a [Circuit](@ref) object. Usage of `MeasurementOperator` can be seen in the [Circuit Example](@ref) shown below.

```@docs
MeasurementOperator
```

## Circuit
[Circuit](@ref) objects combine a `Vector{Moment}` object and `Vector{MeasurementOperator}` object. Applying a `Circuit` to a given input quantum state outputs the various expectation values from the measurement operators defined in the [`Qaintessent.MeasurementOperator`](@ref) objects. A simple circuit is shown in the [Circuit Example](@ref).

```@docs
Circuit{N}
apply(ψ::Vector{<:Complex}, c::Circuit{N}) where {N}
```

### Circuit Example

```jldoctest Circuit
julia> N = 3;

julia> cgs = CircuitGate[circuit_gate(1, HadamardGate()), circuit_gate(2, HadamardGate()), circuit_gate(3, HadamardGate())];

julia> meas = [MeasurementOperator(X, (1,))]; # Measure with regard to X basis on first qubit.

julia> c = Circuit{N}(cgs, meas) # create circuit object

    3 —[H ]—
            
    2 —[H ]—
            
    1 —[H ]—

```
Applying the [Circuit](@ref) object to a 3-qubit quantum state all in the ground 0 state.

```jldoctest Circuit
julia> ψ = ComplexF64[1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];

julia> apply(ψ, c)
1-element Array{Float64,1}:
 0.9999999999999996
```
