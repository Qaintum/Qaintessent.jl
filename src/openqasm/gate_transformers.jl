using Qaintessent

function h(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), HadamardGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), HadamardGate(), N))
    end
end

function x(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), X, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), X, N))
    end
end

function y(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), Y, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), Y, N))
    end
end

function z(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), Z, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), Z, N))
    end
end

function s(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), SGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), SGate(), N))
    end
end

function sdg(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), SdagGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), SdagGate(), N))
    end
end

function t(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), TGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), TGate(), N))
    end
end

function tdg(args, outs...; ccntrl=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), TdagGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), TdagGate(), N))
    end
end

function cx(args, outs...)
    if length(outs) == 2
        qasm_cgc(controlled_circuit_gate((outs[1]), (outs[2]), X, N))
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), X, N))
    end
end

function cz(args, outs...)
    if length(outs) == 2
        qasm_cgc([controlled_circuit_gate((outs[1]),(outs[2]), Z, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), Z, N))
    end
end

function rx(args, outs...)
    if length(outs) == 1
        qasm_cgc(single_qubit_circuit_gate((outs[1]), RxGate(args[1]), N))
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), RxGate(args[1]), N))
    end
end

function ry(args, outs...)
    if length(outs) == 1
        qasm_cgc(single_qubit_circuit_gate((outs[1]), RyGate(args[1]), N))
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), RyGate(args[1]), N))
    end
end

function rz(args, outs...)
    if length(outs) == 1
        qasm_cgc(single_qubit_circuit_gate((outs[1]), RzGate(args[1]), N))
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), RzGate(args[1]), N))
    end
end

function crz(args, outs...)
    if length(outs) == 2
        qasm_cgc(controlled_circuit_gate((outs[1]), (outs[2]), RzGate(args[1]), N))
    else
        qasm_cgc(controlled_circuit_gate((outs[1:end-1]...,), (outs[end]), RzGate(args[1]), N))
    end
end

function u(args, outs...)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), RzGate(args[1]), N),
                single_qubit_circuit_gate((outs[1]), RyGate(args[2]), N),
                single_qubit_circuit_gate((outs[1]), RzGate(args[3]), N)
                ])
    else
        qasm_cgc([controlled_circuit_gate((ccntrl), (outs[1]), RzGate(args[1]), N),
                controlled_circuit_gate(ccntrl, (outs[1]), RyGate(args[2]), N),
                controlled_circuit_gate(ccntrl, (outs[1]), RzGate(args[3]), N)
                ])
    end
end
