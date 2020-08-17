using MLStyle
using Qaintessent


@as_record CircuitGateChain
@as_record Moment
@as_record CircuitGate
@as_record ControlledGate
@as_record HadamardGate
@as_record XGate
@as_record YGate
@as_record ZGate
@as_record RxGate
@as_record RyGate
@as_record RzGate
@as_record SGate
@as_record TGate
@as_record SdagGate
@as_record TdagGate


function trans_circuit(c::Circuit)
    trans_circuit(c.cgc)
end

function trans_bitarray(ba::BitArray{1})
    name_register = "creg creg" * string(classical_register_count) * "[" * string(length(ba)) * "];\n"
    global classical_register_count += 1
    return name_register
end

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
function trans_cg(cg::CircuitGate{M,N,G}) where {M,N,G}
    @match cg begin
        CircuitGate(iwire=iwire,
            gate=gate,
            ccntrl=ccntrl
            ) => let id = trans_g(gate)
                    if !all(isinteger.(ccntrl))
                        nnint = parse(Int64, match(integer_pattern, string(ccntrl[1])).match[1]) - 1
                        return "//if(creg1==" * string(nnint) * ") " * string(id) * join( "qreg" .* string.(iwire) ,",") * "; // If clauses are not yet supported."
                    end
                    return string(id) * " " * join( "qregister[" .* string.(iwire) .* "]" ,",") * ";"
                end
    end
end


function trans_moment(m::Moment{N}) where {N}
    @match m begin
        Moment(
            gates=gates
        ) => let gates=gates
                return reduce(vcat, trans_cg.(gates))
            end
    end
end

function cgc2qasm(cgc::CircuitGateChain{N}) where {N}
    global classical_register_count = 1
    @match cgc begin
        CircuitGateChain(moments=moments,
            creg=creg) => let moments=moments,
                              creg=creg
                                header = """OPENQASM 2.0;\n
                                    include "qelib1.inc";\n
                                    """
                                header *= "qreg qregister[" * string(N) * "];\n"
                                header = header .* trans_bitarray.(creg)
                                header = header .* join(trans_moment.(moments), "\n")
                                return header[1]
                            end
    end
end
