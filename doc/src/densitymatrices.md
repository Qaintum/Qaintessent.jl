# Density Matrix representation

Qaintessent.jl supports the manipulation of mixed states as [DensityMatrix](@ref)(s).

```@meta
CurrentModule = Qaintessent
```

## Density Matrices
```@docs
DensityMatrix
apply(cgs::Vector{<:CircuitGate}, ρ::DensityMatrix)
apply(c::Circuit{N}, ρ::DensityMatrix) where {N}
```

