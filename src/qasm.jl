"""
    custom_gate()

creates template functions from OpenQASM style documentation
"""
function custom_gate!(iwire_dict::Dict, cwire_dict::Dict, gate_dict::Dict, gate_name::String, gate_wires::String, gate_def::String; gate_vars::Union{String, Nothing}=nothing)
    gate_var_expression = r"([a-z][a-zA-Z0-9\_]*)(\(*.*\)*) ([a-z0-9A-Z\[\],]*);"
    if haskey(gate_dict, gate_name)
        return
    end
    gate_wires = strip.(split(gate_wires, ","))

    if !isnothing(gate_vars)
        gate_vars = strip.(split(gate_vars, ","))
    end


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
    gate_dict[gate_name] = gate_func
end

"""
    import_file(filename::String; type="QASM")
converts a .QASM file to a CircuitGateChain
"""
function import_file(filename::String; type="QASM")
    isfile(filename) || error("File '"* filename * "' does not exist")
    h_func =
    function h_func(wires, vars, N)
        if any(wires.<0)
            creg = fill(0, min(abs.(wires)...))
        else
            creg = Int[]
        end
        CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], HadamardGate(), N)]; creg=creg)
    end
    x_func =
    function x_func(wires, vars, N)
        if any(wires.<0)
            creg = fill(0, min(abs.(wires)...))
        else
            creg = Int[]
        end
        CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], XGate(), N)]; creg=creg)
    end
    cx_func =
    function cx_func(wires, vars, N)
        if any(wires.<0)
            creg = fill(0, min(abs.(wires)...))
        else
            creg = Int[]
        end
        CircuitGateChain{N}([controlled_circuit_gate(wires[1], wires[2], XGate(), N)]; creg=creg)
    end
    ccx_func =
    function ccx_func(wires, vars, N)
        if any(wires.<0)
            creg = fill(0, min(abs.(wires)...))
        else
            creg = Int[]
        end
        CircuitGateChain{N}([controlled_circuit_gate((wires[1], wires[2]), wires[3], XGate(), N)]; creg=creg)
    end
    u_func =
    function u_func(wires, vars, N)
        if any(wires.<0)
            creg = fill(0, min(abs.(wires)...))
        else
            creg = Int[]
        end
        CircuitGateChain{N}([
        single_qubit_circuit_gate(wires[1], RxGate(vars[1]), N),
        single_qubit_circuit_gate(wires[1], RyGate(vars[2]), N),
        single_qubit_circuit_gate(wires[1], RzGate(vars[3]), N)
        ]; creg=creg)
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
    gate_novar_expression = r"gate ([a-z][a-zA-z0-9\_]*) (.*)"
    measure_expression = r"measure ([a-z][a-zA-Z0-9\_]*)(?:\[([0-9]+)\])?\s*->\s*([a-z][a-zA-Z0-9\_]*)(?:\[(.*)\])?"
    custom_gate_expression = r"([a-z][a-zA-z0-9\_]*)"
    custom_var_expression = r"([a-z][a-zA-z0-9\_]*)\((.*)\)"
    custom_wire_expression = r"([a-z][a-zA-Z0-9\_]*)(?:\[([0-9]+)\])?"
    open(filename) do f
        Base.mark(f)
        # Initial loop finds size and constructs wire_dict that links qubit registers to iwires in Circuit
        for l in eachline(f)
            line = lstrip(l)
            if startswith(line, "#") || startswith(line, "//")
                continue
            elseif startswith(line, "qreg")
                m = Base.match(register_expression, line)
                name = m.captures[1]
                size = parse(Int64, m.captures[2])
                !haskey(iwire_dict, name) || error("QASM file redefines " * name * " qubit register")
                iwire_dict[name] = collect(N+1:N+size)
                N += size
            elseif startswith(line, "creg")
                m = Base.match(register_expression, line)
                name = m.captures[1]
                size = parse(Int64, m.captures[2])
                !haskey(cwire_dict, name) || error("QASM file redefines " * name * " classical register")
                cwire_dict[name] = collect(M-1:-1:M-size)
                M -= size
            elseif startswith(line, "gate")
                readuntil(f, "{"; keep=false)
                gate_def = readuntil(f, "}"; keep=false)
                m = Base.match(gate_var_expression, line)

                if m == nothing
                    m = Base.match(gate_novar_expression, line)
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

        global cgs = CircuitGateChain{N}(AbstractCircuitGate{N}[]; creg=fill(0, abs(M)))
        global mops = AbstractMatrix[]
        Base.reset(f)
        in_gate_def = false
        for l in eachline(f)
            line = lstrip(l)
            if startswith(line, "#") || startswith(line, "//")
                continue
            elseif startswith(line, "{")
                in_gate_def = true
            elseif startswith(line, "}")
                in_gate_def = false
            elseif startswith(line, "measure")
                m = Base.match(measure_expression, line)
                input = m.captures[1]
                input_wires = m.captures[2]
                output = m.captures[3]
                output_wires = m.captures[4]
                input_register = iwire_dict[input]
                output_register = cwire_dict[output]

                if !isnothing(input_wires)
                    input_wires = [input_register[parse(Int64, input_wires)+1]]
                    input_length = 1
                else
                    input_length = length(input_register)
                end

                if !isnothing(output_wires)
                    output_wires = [output_register[parse(Int64, output_wires)+1]]
                    output_length = 1
                else
                    output_length = length(output_register)
                end

                input_length == output_length || error("Length of registers for measurement must match")
                for i in 1:input_length
                    output_register[i] = input_register[i]
                    m = kron(Matrix{Float64}(I, 2^(input_register[i]-1), 2^(input_register[i]-1)), Qaintessent.matrix(ZGate()))
                    m = kron(m, Matrix{Float64}(I, 2^(N-input_register[i]), 2^(N-input_register[i])))
                    @assert size(m) == (2^N, 2^N)
                    push!(mops, m)
                end
            else
                custom_gate_match = Base.match(custom_gate_expression, line)
                if !isnothing(custom_gate_match) && !in_gate_def
                    custom_gate_name = custom_gate_match.captures[1]
                    register_length = 1
                    if haskey(gate_dict, custom_gate_name)
                        # custom_gate_vars = parse.(Float64, split(custom_gate_match.captures[2], ','))
                        custom_var_match = Base.match(custom_var_expression, line)
                        custom_gate_vars = nothing
                        if !isnothing(custom_var_match) && !isnothing(custom_var_match.captures[2])
                            custom_gate_vars = parse.(Float64, split(custom_var_match.captures[2], ","))
                        end
                        custom_registers = [strip(x, [';', ' ']) for x in split(split(line;limit=2)[2], ",")]
                        custom_wires = []
                        for register in custom_registers
                            m = Base.match(custom_wire_expression, register)
                            name = m.captures[1]
                            reg = m.captures[2]
                            current_register = iwire_dict[name]
                            if !isnothing(reg)
                                push!(custom_wires, iwire_dict[name][parse(Int64,reg)+1])
                            else
                                register_length == 1 || length(iwire_dict[name]) == register_length || error("Register lengths must match when applying gates to multiple wires")
                                register_length = length(iwire_dict[name])
                                push!(custom_wires, iwire_dict[name])
                            end
                        end
                        for i in 1:register_length
                            custom_gate_wires = []
                            for j in 1:length(custom_wires)
                                if length(custom_wires[j]) > 1
                                    push!(custom_gate_wires, custom_wires[j][i])
                                else
                                    push!(custom_gate_wires, custom_wires[j][1])
                                end
                            end
                            cgs *= gate_dict[custom_gate_name](custom_gate_wires, custom_gate_vars, N)
                        end
                    end
                end

            end
        end
    end
    return Circuit{N}(cgs, MeasurementOps{N}(mops))
end

"""
    export_file(filename::String; type="QASM")
CircuitGateChain or Circuit to a QASM file
"""
function export_file(circuit::Circuit{N}; type="QASM") where {N}

end
