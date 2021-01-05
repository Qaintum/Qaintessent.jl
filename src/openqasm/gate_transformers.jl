
function custom_gate_h(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), HadamardGate())])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), HadamardGate(), (outs[1:end-1]...,)))
    end
end

function custom_gate_x(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), X)])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), X, (outs[1:end-1]...,)))
    end
end

function custom_gate_y(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), Y)])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), Y, (outs[1:end-1]...,)))
    end
end

function custom_gate_z(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), Z)])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), Z, (outs[1:end-1]...,)))
    end
end

function custom_gate_s(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), SGate())])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), SGate(), (outs[1:end-1]...,)))
    end
end

function custom_gate_sdg(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), SdagGate())])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), SdagGate(), (outs[1:end-1]...,)))
    end
end

function custom_gate_t(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), TGate())])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), TGate(), (outs[1:end-1]...,)))
    end
end

function custom_gate_tdg(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate((outs[1]), TdagGate())])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), TdagGate(), (outs[1:end-1]...,)))
    end
end

function custom_gate_cx(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        append!(qasm_cgc, circuit_gate((outs[2]), X, (outs[1])))
    else
        append!(qasm_cgc, circuit_gate((outs[end]), X, (outs[1:end-1]...,)))
    end
end

function custom_gate_cz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        append!(qasm_cgc, [circuit_gate((outs[2]), Z, (outs[1]))])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), Z, (outs[1:end-1]...,)))
    end
end

function custom_gate_rx(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, circuit_gate((outs[1]), RxGate(args[1])))
    else
        append!(qasm_cgc, circuit_gate((outs[end]), RxGate(args[1]), (outs[1:end-1]...,)))
    end
end

function custom_gate_ry(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, circuit_gate((outs[1]), RyGate(args[1])))
    else
        append!(qasm_cgc, circuit_gate((outs[end]), RyGate(args[1]), (outs[1:end-1]...,)))
    end
end

function custom_gate_rz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, circuit_gate((outs[1]), RzGate(args[1])))
    else
        append!(qasm_cgc, circuit_gate((outs[end]), RzGate(args[1]), (outs[1:end-1]...,)))
    end
end

function custom_gate_crz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        append!(qasm_cgc, circuit_gate((outs[2]), RzGate(args[1]), (outs[1])))
    else
        append!(qasm_cgc, circuit_gate((outs[end]), RzGate(args[1]), (outs[1:end-1]...,)))
    end
end

function custom_gate_u(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 1
        append!(qasm_cgc, [circuit_gate(outs[1], RzGate(args[1])),
                  circuit_gate(outs[1], RyGate(args[2])),
                  circuit_gate(outs[1], RzGate(args[3]))
                ])
    else
        append!(qasm_cgc, [circuit_gate(outs[1], RzGate(args[1], ccntrl)),
                  circuit_gate(outs[1], RyGate(args[2]), ccntrl),
                  circuit_gate(outs[1], RzGate(args[3]), ccntrl)
                ])
    end
end
