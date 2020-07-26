# OpenQASM

Qaintessent.jl provides some basic OpenQASM 2.0 support. This allows the importing of most QASM files as Circuit objects. Currently, only basic exports of Circuits is supported: custom measurements and multi-control CircuitGates cannot be exported. An example of OpenQASM usage is seen in the example below.

```@meta
CurrentModule = Qaintessent
```

```@docs
    import_qasm(filename::String)
```

```@docs
    import_file(filename::String; type="QASM")
```

```@docs
    export_qasm(circuit::Circuit{N}, filename::String) where {N}
```

```@docs
    export_file(circuit::Circuit{N}, filename::String; type="QASM") where {N}
```

#### Example
```jldoctest
using Qaintessent 

N = 3
filename = "reference_circuit.qasm"

cgc_ref = CircuitGateChain{N}([
    single_qubit_circuit_gate(3, XGate(), N),
    two_qubit_circuit_gate(3, 2, SwapGate(), N),
    single_qubit_circuit_gate(1, RyGate(0.3π), N),
    single_qubit_circuit_gate(2, RzGate(2.4), N),
    single_qubit_circuit_gate(2, TGate(), N),
    controlled_circuit_gate(1,2, XGate(), N),
    controlled_circuit_gate(3,1, RzGate(0.3), N)
    ])

meas = MeasurementOps{N}(AbstractMatrix[])
c_ref = Circuit{N}(cgc_ref, meas)
export_file(c_ref, filename)

c = import_file(filename)
rm(filename)

c
# output


    1 ————————————————————•———
                          |
    2 ————————x————[T ]——[X ]—
              |
    3 —[X ]———x———————————————
```