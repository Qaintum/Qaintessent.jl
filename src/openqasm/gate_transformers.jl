
function custom_gate_cx(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        append!(qasm_cgc, [circuit_gate((outs[2]), X, (outs[1]))])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), X, (outs[1:end-1]...,)))
    end
end

function custom_gate_cy(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        append!(qasm_cgc, [circuit_gate((outs[2]), Y, (outs[1]))])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), Y, (outs[1:end-1]...,)))
    end
end

function custom_gate_cz(args, outs...; qasm_cgc=nothing, N=nothing)
    if length(outs) == 2
        append!(qasm_cgc, [circuit_gate((outs[2]), Z, (outs[1]))])
    else
        append!(qasm_cgc, circuit_gate((outs[end]), Z, (outs[1:end-1]...,)))
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
