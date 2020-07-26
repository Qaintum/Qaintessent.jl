using Formatting

"""
    custom_gate()

creates template functions from OpenQASM style documentation
"""
function custom_gate!(iwire_dict::Dict, cwire_dict::Dict, gate_dict::Dict, gate_name::String, gate_wires::String, gate_def::String; gate_vars::Union{String, Nothing}=nothing)
    gate_var_expression = r"([a-zA-Z][a-zA-Z0-9\_]*)(\(*.*\)*) ([a-z0-9A-Z\[\],]*);"
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
                            replace(local_var, "pi"=>"π")
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
converts an external file to a Circuit object
"""
function import_file(filename::String; type="QASM")
    if type == "QASM"
        return import_qasm(filename::String)
    end
end

h_func =
function h_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], HadamardGate(), N)])
end
t_func =
function t_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], TGate(), N)])
end
tdag_func =
function tdag_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], TdagGate(), N)])
end
s_func =
function s_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], SGate(), N)])
end

sdag_func =
function sdag_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], SdagGate(), N)])
end

x_func =
function x_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], XGate(), N)])
end

y_func =
function y_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], YGate(), N)])
end

z_func =
function z_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], ZGate(), N)])
end

rx_func =
function rx_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], RxGate(vars[1]), N)])
end

ry_func =
function ry_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], RyGate(vars[1]), N)])
end

rz_func =
function rz_func(wires, vars, N)
    CircuitGateChain{N}([single_qubit_circuit_gate(wires[1], RzGate(vars[1]), N)])
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
    CircuitGateChain{N}([
    single_qubit_circuit_gate(wires[1], RzGate(vars[2]), N),
    single_qubit_circuit_gate(wires[1], RyGate(vars[1]), N),
    single_qubit_circuit_gate(wires[1], RzGate(vars[3]), N)
    ])
end

cu1_func =
function cu1_func(wires, vars, N)
    CircuitGateChain{N}([
    single_qubit_circuit_gate(wires[1], RzGate(vars[1]/2), N),
    controlled_circuit_gate(wires[1], wires[2], XGate(), N),
    single_qubit_circuit_gate(wires[1], RzGate(-vars[1]/2), N),
    controlled_circuit_gate(wires[1], wires[2], XGate(), N),
    single_qubit_circuit_gate(wires[1], RzGate(vars[1]/2), N),
    ])
end

swap_func =
function swap_func(wires, vars, N)
    CircuitGateChain{N}([two_qubit_circuit_gate(wires[1], wires[2], SwapGate(), N)])
end

ref_gate_dict = Dict("x"=>x_func,
                 "y"=>y_func,
                 "z"=>z_func,
                 "s"=>s_func,
                 "sdag"=>sdag_func,
                 "t"=>t_func,
                 "tdag"=>tdag_func,
                 "h"=>h_func,
                 "swap"=>swap_func,
                 "cx"=>cx_func,
                 "CX"=>cx_func,
                 "ccx"=>ccx_func,
                 "CCX"=>ccx_func,
                 "U"=>u_func,
                 "cu1"=>cu1_func
                 )

"""
    import_file(filename::String; type="QASM")
converts a .QASM file to a CircuitGateChain
"""
function import_qasm(filename::String)
    gate_dict = deepcopy(ref_gate_dict)
    isfile(filename) || error("File '"* filename * "' does not exist")
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
            elseif startswith(line, "if")
                error("Qaintessent.jl currently does not support the `if (creg===x)` command in OpenQASM`")
            elseif startswith(line, "reset")
                error("Qaintessent.jl currently does not support the `reset` command in OpenQASM`")
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
                    input_wires = input_register[parse(Int64, input_wires)+1]
                    input_length = 1
                else
                    input_wires = input_register[1] # In case register is size 1
                    input_length = length(input_register)
                end

                if !isnothing(output_wires)
                    output_wires = parse(Int64, output_wires)+1
                    output_length = 1
                else
                    output_wires = 1 # In case register is size 1
                    output_length = length(output_register)
                end

                input_length == output_length || error("Length of registers for measurement must match")

                if output_length > 1
                    for i in 1:input_length
                        output_register[i] = input_register[i]
                        m = kron(Matrix{ComplexF64}(I, 2^(input_register[i]-1), 2^(input_register[i]-1)), Qaintessent.matrix(ZGate()))
                        m = kron(m, Matrix{ComplexF64}(I, 2^(N-input_register[i]), 2^(N-input_register[i])))

                        @assert size(m) == (2^N, 2^N)
                        push!(mops, m)
                    end
                else
                    output_register[output_wires] = input_wires
                    m = kron(Matrix{ComplexF64}(I, 2^(input_wires-1), 2^(input_wires-1)), Qaintessent.matrix(ZGate()))
                    m = kron(m, Matrix{ComplexF64}(I, 2^(N-input_wires), 2^(N-input_wires)))
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
                            custom_gate_vars = replace(custom_var_match.captures[2], "pi"=>"π")
                            custom_gate_vars = eval.(Meta.parse.((split(custom_gate_vars, ","))))
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

gate_name_dict = Dict(
                    XGate => "x {1};\n",
                    YGate => "y {1};\n",
                    ZGate => "z {1};\n",
                    RxGate => "rx({2}) {1};\n",
                    RyGate => "ry({2}) {1};\n",
                    RzGate => "rz({2}) {1};\n",
                    TGate => "t {1};\n",
                    TdagGate => "tdag {1};\n",
                    SGate => "s {1};\n",
                    SdagGate => "sdag {1};\n",
                    HadamardGate => "h {1};\n",
                    SwapGate => "swap {1};\n",
                        )

controlled_gate_name_dict = Dict(
                            XGate => "cx {1};\n",
                            YGate => "cy {1};\n",
                            ZGate => "cz {1};\n",
                            HadamardGate => "ch {1};\n",
                            RzGate => "crz({2}) {1};\n",
                            )
"""
    export_qasm(filename::String; type="QASM")
Circuit to a QASM file
"""
function export_qasm(circuit::Circuit{N}, filename::String) where {N}
    output = """OPENQASM 2.0;\ninclude "qelib1.inc";\n\n"""
    M = length(circuit.cgc.creg)
    if M > 0
        output *= format("//Defining registers\nqreg q[{1}];\ncreg c[{2}];\n\n", N, M)
    else
        output *= format("//Defining registers\nqreg q[{1}];\n\n", N)
    end

    for moment in circuit.cgc
        for circuit_gate in moment
            if isa(circuit_gate.gate, ControlledGate)
                if length(get_controls(circuit_gate)[1]) > 1
                    error(format("Gate {1} has {2} control wires. Multi-control wires are currently not supported for the OpenQASM format.", circuit_gate, length(get_controls(circuit_gate)[1])))
                end
                abstract_gate = circuit_gate.gate.U
                gate_type = typeof(abstract_gate)
                fields = fieldnames(gate_type)

                if haskey(controlled_gate_name_dict, gate_type)
                    frmt = controlled_gate_name_dict[gate_type]
                else
                    error("Controlled Gate " * string(gate_type) * " conversion to OpenQASM is currently not supported in Qaintessent.jl")
                end
            else
                abstract_gate = circuit_gate.gate
                gate_type = typeof(abstract_gate)
                fields = fieldnames(gate_type)
                if haskey(gate_name_dict, gate_type)
                    frmt = gate_name_dict[gate_type]
                else
                    error("Gate " * string(gate_type) * " convertsion to OpenQASM is currently not supported in Qaintessent.jl")
                end
            end

            vars = ""
            for field in fields
                vars *= string(getfield(abstract_gate, field)[1]) * ","
            end
            vars = vars[1:end-1]

            wires = format("q[{1}]", circuit_gate.iwire[1]-1)
            for wire in circuit_gate.iwire[2:end]
                wires *= format(",q[{1}]", wire-1)
            end

            output *= format(frmt, wires, vars)

            end
        end
    open(filename, "w+") do f
        write(f, output)
    end
end

"""
    export_file(filename::String; type="QASM")
Circuit to a export file
"""
function export_file(circuit::Circuit{N}, filename::String; type="QASM") where {N}
    if type == "QASM"
        export_qasm(circuit, filename)
    end
end
