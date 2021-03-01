# Circuit Gates

A [CircuitGate](@ref) object consists of an [`Qaintessent.AbstractGate`](@ref) and wire ordering information. This allows an arbitrary quantum operation can be encoded as vector of [CircuitGate](@ref) objects and applied to a given quantum state.

```@meta
CurrentModule = Qaintessent
DocTestSetup = quote
    using Qaintessent
end
```

## CircuitGate
```@docs
CircuitGate
circuit_gate
apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}
apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})
sparse_matrix(cg::CircuitGate{M,G}, N::Integer=0) where {M,G <: AbstractGate}
sparse_matrix(cgs::Vector{<:CircuitGate}, N::Integer=0)
```

## Compilation of abitrary unity operation
An arbitrary unitary operation can be compiled into a vector of [CircuitGate](@ref) objects

```@docs
unitary2circuit
```