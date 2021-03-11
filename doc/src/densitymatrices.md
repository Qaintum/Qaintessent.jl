# Density Matrix representation

Qaintessent.jl supports the manipulation of mixed states as [Density Matrices](@ref)(s).

```@meta
CurrentModule = Qaintessent
```

## Density Matrices
Density matrices are stored as linear combination of Pauli matrices. Given a quantum state of `N` qubits, we require a `Vector{<:Real}` of length 4^N to 
store the mixed state.

```@docs
DensityMatrix(v::AbstractVector{<:Real}, N::Integer)
DensityMatrix(v::AbstractVector{<:Real})
```
Several helper functions are added to crate [`DensityMatrix`](@ref) objects from quantum states in statevector form.

```@docs
density_from_statevector(ψ::Vector)
density_from_matrix(ρmat::AbstractMatrix)
```
Once created, they can be treated similar to a statevector, see [Circuit Construction and Usage](@ref).

```@docs
apply(ρ::DensityMatrix, cgs::Vector{<:CircuitGate})
apply(ρ::DensityMatrix, c::Circuit{N}) where {N}
```