using MLStyle


@as_record Moment
@as_record CircuitGate
@as_record ControlledGate
@as_record HadamardGate
@as_record XGate
@as_record YGate
@as_record ZGate
@as_record SwapGate
@as_record RxGate
@as_record RyGate
@as_record RzGate
@as_record SGate
@as_record TGate
@as_record SdagGate
@as_record TdagGate

# function trans_bitarray(ba::BitArray{1})
#     name_register = "creg creg" * string(classical_register_count) * "[" * string(length(ba)) * "];\n"
#     global classical_register_count += 1
#     return name_register
# end

function trans_g(g::AbstractGate)
    @match g begin
        HadamardGate() => "h"
        XGate() => "x"
        YGate() => "y"
        ZGate() => "z"
        SGate() => "s"
        SdagGate() => "sdg"
        TGate() => "t"
        TdagGate() => "tdg"
        SwapGate() => "swap"
        RxGate(θ=θ) => "rx(" * string(θ[1]) * ")"
        RyGate(θ=θ) => "ry(" * string(θ[1]) * ")"
        RzGate(θ=θ) => "rz(" * string(θ[1]) * ")"
        ControlledGate(U=U
            ) => let U=trans_g(U)
                    return "c" * U
                end
    end
end

integer_pattern = r"[0-9]+$"
function trans_cg(cg::CircuitGate{M,G}) where {M,G}
    @match cg begin
        CircuitGate(iwire=iwire,
            gate=gate
            ) => let id = trans_g(gate)
                    return string(id) * " " * join( "qregister[" .* string.(reverse(iwire.-1)) .* "]" ,",") * ";"
                end
    end
end


function trans_moment(m::Moment)
    @match m begin
        Moment(
            gates=gates
        ) => let gates=gates
                expr = join(trans_cg.(gates), "\n")
                return expr
            end
    end
end

"""
    cgc2qasm(c::Circuit{N}) where {N}

converts Circuit{N} object to OpenQASM 2.0 representation
"""
function cgc2qasm(c::Circuit{N}) where {N}
    global classical_register_count = 1
    header = """OPENQASM 2.0;\n
    include "qelib1.inc";\n
    """
    header *= "qreg qregister[" * string(N) * "];\n"
    header = header .* join(trans_moment.(c.moments), "\n")    
    return header
end
