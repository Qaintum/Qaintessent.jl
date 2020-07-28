# Qaintessent.jl Documentation

## Table of Contents
```@contents
Pages = ["index.md", "gates.md", "circuit.md", "gradients.md", "models.md", "view.md", "qasm.md"]
```

```@meta
CurrentModule = Qaintessent
```

## Features

Qaintessent.jl is a digital quantum circuit toolbox and simulator, using Julia's type system to represent quantum circuits symbolically. Quantum circuits are created from combinations of basic quantum gates (See [Gates](@ref) for details regarding elementary gates and [Circuit Construction and Usage](@ref) for circuit construction). Qaintessent.jl represents quantum states in state vector form, or (alternatively) density matrices with respect to the Pauli basis. It also supports gradient calculation for parametrized gates to be used for optimization algorithms, such as Machine Learning or QAOA. Integration with [Flux](https://fluxml.ai) is provided by [Qaintellect.jl](https://github.com/Qaintum/Qaintellect.jl).

## Index

```@index
```
