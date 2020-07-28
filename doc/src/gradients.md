# Gradients

Qaintessent.jl can perform a backward pass through a circuit to compute gradients of a (fictitious) cost function with respect to the circuit parameters and input wavefunction. Since quantum gates are unitary, the intermediate quantum states can be re-computed on the fly during backwards traversal, i.e., need not be stored during the forward pass.

See also [this article](https://fluxml.ai/Zygote.jl/latest/complex) regarding complex differentiation.

```@meta
CurrentModule = Qaintessent
```

```@docs
Qaintessent.gradients
```
