# Gates

Any quantum operation can be represented as a unitary matrix, which are referred to as `gates`. In Qaintessent.jl, basic gates are implemented using the abstract base struct `AbstractGate`. The available basic quantum gates are seen below. Note that wire index information is still needed to when applying such gates. As such, these gates are permuted when creating the basic building blocks of Qaintessent.jl Quantum [Circuit](@ref)(s): [CircuitGate](@ref)(s).

```@meta
CurrentModule = Qaintessent
```
```@docs
AbstractGate
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
## Parametrized Single Qubit Gates
```@docs
RxGate
RyGate
RzGate
RotationGate
PhaseShiftGate
```

## Two Qubit Gates
```@docs
SwapGate
```

## Quantum Gate Matrix Representation
The `Qaintessent.matrix` function can be used to convert `AbstractGate` objects to their matrix representation. 

```@example
using Qaintessent

x = XGate()
Qaintessent.matrix(x)
```


