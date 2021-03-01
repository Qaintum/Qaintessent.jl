# Density Matrix representation

Qaintessent.jl supports the manipulation of mixed states as [DensityMatrix](@ref)(s).

```@meta
CurrentModule = Qaintessent
```

## Density Matrices
```@docs
DensityMatrix
apply(cg::CircuitGate{M,G}, ρ::DensityMatrix) where {M,G}
apply(cg::Vector{<:CircuitGate}, ρ::DensityMatrix)
apply(c::Circuit{N}, ρ::DensityMatrix) where {N}
```

