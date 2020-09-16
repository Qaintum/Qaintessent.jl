
function custom_gate_h(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), HadamardGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), HadamardGate(), N[]))
    end
end

function custom_gate_x(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), X, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), X, N[]))
    end
end

function custom_gate_y(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), Y, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), Y, N[]))
    end
end

function custom_gate_z(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), Z, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), Z, N[]))
    end
end

function custom_gate_s(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), SGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), SGate(), N[]))
    end
end

function custom_gate_sdg(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), SdagGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), SdagGate(), N[]))
    end
end

function custom_gate_t(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), TGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), TGate(), N[]))
    end
end

function custom_gate_tdg(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate((outs[1]), TdagGate(), N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), TdagGate(), N[]))
    end
end

function custom_gate_cx(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        qasm_cgc(controlled_circuit_gate((outs[2]), (outs[1]), X, N[]))
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), X, N[]))
    end
end

function custom_gate_cz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        qasm_cgc([controlled_circuit_gate((outs[2]), (outs[1]), Z, N)])
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), Z, N[]))
    end
end

function custom_gate_rx(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc(single_qubit_circuit_gate((outs[1]), RxGate(args[1]), N[]))
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), RxGate(args[1]), N[]))
    end
end

function custom_gate_ry(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc(single_qubit_circuit_gate((outs[1]), RyGate(args[1]), N[]))
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), RyGate(args[1]), N[]))
    end
end

function custom_gate_rz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc(single_qubit_circuit_gate((outs[1]), RzGate(args[1]), N[]))
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), RzGate(args[1]), N[]))
    end
end

function custom_gate_crz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        qasm_cgc(controlled_circuit_gate((outs[2]), (outs[1]), RzGate(args[1]), N[]))
    else
        qasm_cgc(controlled_circuit_gate((outs[end]), (outs[1:end-1]...,), RzGate(args[1]), N[]))
    end
end

function custom_gate_u(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        qasm_cgc([single_qubit_circuit_gate(outs[1], RzGate(args[1]), N),
                  single_qubit_circuit_gate(outs[1], RyGate(args[2]), N),
                  single_qubit_circuit_gate(outs[1], RzGate(args[3]), N)
                ])
    else
        qasm_cgc([controlled_circuit_gate(outs[1], ccntrl, RzGate(args[1]), N),
                  controlled_circuit_gate(outs[1], ccntrl, RyGate(args[2]), N),
                  controlled_circuit_gate(outs[1], ccntrl, RzGate(args[3]), N)
                ])
    end
end
