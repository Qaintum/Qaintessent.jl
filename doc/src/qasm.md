# OpenQASM

Qaintessent.jl provides some basic OpenQASM 2.0 support. This allows the importing of most QASM files as Circuit objects. Currently, only basic exports of Circuits is supported: custom measurements and multi-control CircuitGates cannot be exported. An example of OpenQASM usage is seen in the example below.

```@meta
CurrentModule = Qaintessent
```

```@docs
    qasm2cgc(txt::String)
```

```@docs
    cgc2qasm(c::Circuit{N}) where {N}
```

#### Example
```jldoctest
using Qaintessent
using LinearAlgebra

N = 3
cgs_ref = CircuitGate[
    circuit_gate(3, XGate()),
    circuit_gate(1, RyGate(0.3π)),
    circuit_gate(2, RzGate(2.4)),
    circuit_gate(2, TGate()),
    circuit_gate(1, XGate(), 2),
    circuit_gate(3, RzGate(0.3), 1)
    ]

meas = MeasurementOperator(Matrix{Float64}(I, 2^N, 2^N), (1,2,3))
c_ref = Circuit{N}(cgs_ref, [meas])
qasm_rep = cgc2qasm(c_ref)

cgc_from_qasm = qasm2cgc(qasm_rep)

cgc_from_qasm
# output
    
    3 —[X ]——————————————[Rz]—
                          |
    2 —[Rz]——[T ]———•—————————
                    |     |
    1 —[Ry]————————[X ]———•———
```
