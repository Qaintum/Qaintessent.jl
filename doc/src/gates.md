# Gates

Any quantum operation can be represented as a unitary matrix, which are referred to as `gates`. In Qaintessent.jl, basic gates are implemented using the abstract struct `AbstractGate`. The available basic quantum gates are seen below. Note that at this level of abstraction, wire index information is ignored and matrixes are permuted such that the fastest running qubit is also the fastest running index. 

These gates are permuted when creating the basic building blocks of Qaintessent.jl Quantum [Circuit](@ref)(s): [CircuitGate](@ref)(s).

```@meta
CurrentModule = Qaintessent
```
```@docs
AbstractGate
```

## AbstractGate helper functions

```@docs
 matrix(::AbstractGate)
 sparse_matrix(::AbstractGate)
```

```@docs
Base.adjoint(::AbstractGate)
num_wires(::AbstractGate)
Base.isapprox(g1::G, g2::G) where {G <: AbstractGate}
```

```@docs
kron(m::Union{UniformScaling{Bool},G}...) where {G<:AbstractGate}
```

## Basic Single Qubit Gates

```@docs
XGate
YGate
ZGate
HadamardGate
SGate
SdagGate
TGate
TdagGate
```

## Basic Two Qubit Gates
```@docs
SwapGate

```

## Parametrized Gates
```@docs
RxGate
RyGate
RzGate
RotationGate
PhaseShiftGate
EntanglementXXGate
EntanglementYYGate
EntanglementZZGate
```

## Controlled Gates
Arbitrary controlled gates using the basic gates defined above are available.

```@docs
ControlledGate
```

## ControlledGate helper functions

```@docs
control_wires
target_wires
```