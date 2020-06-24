"""
    custom_gate()

creates template functions from OpenQASM style documentation
"""
function custom_gate!(iwire_dict::Dict, cwire_dict::Dict, gate_dict::Dict, gate_name::String, gate_wires::String, gate_def::String; gate_vars::Union{String, Nothing}=nothing)
    gate_var_expression = r"([a-z][a-zA-Z0-9\_]*)(\(*.*\)*) ([a-z\[\],]*);"
    if haskey(gate_dict, gate_name)
        return
    end
    gate_wires = strip.(split(gate_wires, ","))

    if !isnothing(gate_vars)
        gate_vars = strip.(split(gate_vars, ","))
    end
    # println(gate_def)
    # println(gate_vars)

    function gate_func(wires, vars, N)
        cgs = CircuitGateChain{N}(AbstractCircuitGate{N}[])
        local_iwire_dict = Dict()
        local_vars_dict = Dict()
        for (i, wire) in enumerate(gate_wires)
            local_iwire_dict[wire] = wires[i]
        end

        if !isnothing(gate_vars)
            for (j, var) in enumerate(gate_vars)
                local_vars_dict[var] = vars[j]
            end
        end
        for line in split(gate_def, "\n")

            m = Base.match(gate_var_expression, line)
            if m != nothing
                name = m.captures[1]
                rvars = strip.(split(m.captures[2][2:end-1], ","))
                if rvars == [""]
                    rvars = nothing
                end
                rwires = strip.(split(m.captures[3], ","))
                local_wires = []
                local_vars = []
                for local_wire in rwires
                    push!(local_wires, local_iwire_dict[local_wire])
                end

                if !isnothing(gate_vars) && !isnothing(rvars)
                    for local_var in rvars
                        if haskey(local_vars_dict, local_var)
                            push!(local_vars, local_vars_dict[local_var])
                        else
                            push!(local_vars, parse(Float64, local_var))
                        end
                    end
                end

                haskey(gate_dict, name) || error(name * " has not been defined!")
                cgs *= gate_dict[name](local_wires, local_vars, N)
            end
        end
        return cgs
    end
    wires = [1,2,3]
    vars = [0.2, π]
    c = gate_func(wires, vars, 4)
    println(c)
    gate_dict[gate_name] = gate_func
end

"""
    import_file(filename::String; type="QASM")
converts a .QASM file to a CircuitGateChain
"""
function import_file(filename::String; type="QASM")
    isfile(filename) || error("File '"* filename * "' does not exist")
    h_func =
    function h_func(wires, vars, N; ccntrl=Int[])
        CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], HadamardGate(), N; ccntrl=ccntrl)])
    end
    x_func =
    function x_func(wires, vars, N; ccntrl=Int[])
        CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], XGate(), N; ccntrl=ccntrl)])
    end
    cx_func =
    function cx_func(wires, vars, N; ccntrl=Int[])
        CircuitGateChain{N}([controlled_circuit_gate(wires[1], wires[2], XGate(), N; ccntrl=ccntrl)])
    end
    ccx_func =
    function ccx_func(wires, vars, N; ccntrl=Int[])
        CircuitGateChain{N}([controlled_circuit_gate((wires[1], wires[2]), wires[3], XGate(), N; ccntrl=ccntrl)])
    end
    u_func =
    function u_func(wires, vars, N; ccntrl=Int[])
        CircuitGateChain{N}([
        single_qubit_circuit_gate(wires[1], RxGate(vars[1]), N; ccntrl=ccntrl),
        single_qubit_circuit_gate(wires[1], RyGate(vars[2]), N; ccntrl=ccntrl),
        single_qubit_circuit_gate(wires[1], RyGate(vars[3]), N; ccntrl=ccntrl)
        ])
    end
    gate_dict = Dict("x"=>x_func,
                     "X"=>x_func,
                     "h"=>h_func,
                     "H"=>h_func,
                     "cx"=>cx_func,
                     "CX"=>cx_func,
                     "ccx"=>ccx_func,
                     "CCX"=>ccx_func,
                     "u"=>u_func,
                     "U"=>u_func,
                     )
    iwire_dict = Dict()
    cwire_dict = Dict()
    N = 0
    M = 0
    creg = []
    register_expression = r"([a-z\_][a-zA-z0-9\_]*)\[([0-9]+)\]"
    gate_var_expression = r"gate ([a-z][a-zA-z0-9\_]*)\((.*)\) (.*)"
    gate__novar_expression = r"gate ([a-z][a-zA-z0-9\_]*) (.*)"
    open(filename) do f
        # Initial loop finds size and constructs wire_dict that links qubit registers to iwires in Circuit
        for l in eachline(f)
            line = lstrip(l)
            if startswith(line, "#")
                continue
            elseif startswith(line, "qreg")
                m = Base.match(register_expression, line)
                name = m.captures[1]
                size = parse(Int64, m.captures[2])
                !haskey(iwire_dict, name) || error("QASM file redefines " * name * " qubit register")
                iwire_dict[name] = [N+1:N+size]
                N += size
            elseif startswith(line, "creg")
                m = Base.match(register_expression, line)
                name = m.captures[1]
                size = parse(Int64, m.captures[2])
                !haskey(cwire_dict, name) || error("QASM file redefines " * name * " classical register")
                cwire_dict[name] = [M-1:M-size]
                M -= size
            elseif startswith(line, "gate")
                readuntil(f, "{"; keep=false)
                gate_def = readuntil(f, "}"; keep=false)
                m = Base.match(gate_var_expression, line)
                if m == nothing
                    m = Base.match(gate__novar_expression, line)
                    gate_name = String(m[1])
                    gate_vars = nothing
                    gate_wires = String(m[2])
                else
                    gate_name = String(m[1])
                    gate_vars = String(m[2])
                    gate_wires = String(m[3])
                end
                custom_gate!(iwire_dict, cwire_dict, gate_dict, gate_name, gate_wires, gate_def; gate_vars=gate_vars)
            end
        end
        for l in eachline(f)
            line = lstrip(l)
            if startswith(line, "#")
                continue
            end
        end
    end
end

"""
    export_file(filename::String; type="QASM")
CircuitGateChain or Circuit to a QASM file
"""
function export_file(circuit::Circuit{N}; type="QASM") where {N}

end
