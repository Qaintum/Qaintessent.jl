using RBNF
using MLStyle
using PrettyPrint


@as_record Token
@as_record Struct_mainprogram
@as_record Struct_ifstmt
@as_record Struct_gate
@as_record Struct_gatedecl
@as_record Struct_decl
@as_record Struct_barrier_ids
@as_record Struct_reset
@as_record Struct_measure
@as_record Struct_iduop
@as_record Struct_u
@as_record Struct_cx
@as_record Struct_h
@as_record Struct_x
@as_record Struct_y
@as_record Struct_z
@as_record Struct_s
@as_record Struct_sdg
@as_record Struct_t
@as_record Struct_tdg
@as_record Struct_rx
@as_record Struct_ry
@as_record Struct_rz
@as_record Struct_ch
@as_record Struct_crz
@as_record Struct_idlist
@as_record Struct_mixedlist
@as_record Struct_argument
@as_record Struct_explist
@as_record Struct_pi
@as_record Struct_fnexp
@as_record Struct_neg
@as_record Struct_bin

function rec(ctx_tokens)
    function app(op, args...)
        args = map(rec, args)
        op = Symbol(op)
        :($op($(args...)))
    end

    @match ctx_tokens begin
        Struct_pi(_) => Base.pi
        Token{:id}(str=str) => Symbol(str)
        Token{:nnreal}(str=str) => parse(Float64, str)
        Token{:nninteger}(str=str) => parse(Int64, str)
        Struct_neg(value=value) => :(-$(rec(value)))
        Struct_bin(l=l, op=Token(str=op), r=r) => app(op, l, r)
        Struct_idlist(hd=Token(str=hd), tl=nothing) => [Symbol(hd)]
        Struct_idlist(hd=Token(str=hd), tl=tl) => [Symbol(hd), rec(tl)...]

        Struct_explist(hd=hd, tl=nothing) => [rec(hd)]
        Struct_explist(hd=hd, tl=tl) => [rec(hd), rec(tl)...]

        Struct_mixedlist(hd=hd, tl=nothing) => [rec(hd)]
        Struct_mixedlist(hd=hd, tl=tl) => [rec(hd), rec(tl)...]

        Struct_argument(id=Token(str=id), arg=nothing) => Symbol(id)
        Struct_argument(id=Token(str=id), arg=Token(str=int)) =>
            let ind = parse(Int, int) + 1 # due to julia 1-based index
                :($(Symbol(id))[$ind])
            end
        Struct_fnexp(fn = Token(str=fn), arg=arg) =>
            let fn = @match fn begin
                        "sin" => sin
                        "cos" => cos
                        "tan" => tan
                        "exp" => exp
                        "ln"  => log
                        "sqrt"=> sqrt
                        _     => error("not impl yet")
                    end
                app(fn, arg)
            end
        _ => nothing
    end
end

function trans_reg(ctx_tokens, cregs, qregs)
    function app(op, args...)
        args = map(rec, args)
        op = Symbol(op)
        :($op($(args...)))
    end

    @match ctx_tokens begin
        Struct_decl(
            regtype = Token(str=regtype),
            id = Token(str=id),
            int = Token(str = n)
        ) =>
            let id = Symbol(id),
                n = parse(Int, n)
                if regtype == "qreg"
                    return :($id = qreg($n); push!($qregs, $id))
                else
                    return :($id = creg($n); push!($cregs, $id))
                end
            end

        Struct_mainprogram(
            prog = stmts
        ) =>
            let stmts = trans_reg.(stmts, cregs, qregs)
                stmts
            end
        _ => nothing
    end
end

function trans_gates(ctx_tokens, qasm_cgc,  N)
    function app(op, args...)
        args = map(rec, args)
        op = Symbol(op)
        :($op($(args...)))
    end
    @match ctx_tokens begin
        Struct_iduop(gate_name = Token(str=gate_name), args=nothing, outs=outs) =>
            let refs = rec(outs),
                gate_name = Symbol("custom_gate_"*gate_name)
                :($gate_name(nothing, $(refs...); qasm_cgc=$qasm_cgc, N=$N))
            end

        Struct_iduop(gate_name = Token(str=gate_name), args=exprlist, outs=outs) =>
            let refs = rec(outs),
                exprs = Expr(:tuple, rec(exprlist)...),
                gate_name = Symbol("custom_gate_"*gate_name)
                :($gate_name($exprs, $(refs...); qasm_cgc=$qasm_cgc, N=$N))
            end

        Struct_cx(out1=out1, out2=out2) =>
            let ref1 = rec(out1),
                ref2 = rec(out2)
                :($qasm_cgc([controlled_circuit_gate(($ref2), ($ref1), X, $N[])]))
            end

        Struct_ch(out1=out1, out2=out2) =>
            let ref1 = rec(out1),
                ref2 = rec(out2)
                :($qasm_cgc([controlled_circuit_gate(($ref2), ($ref1), HadamardGate(), $N[])]))
            end

        Struct_u(in1=in1, in2=in2, in3=in3, out=out) =>
            let (a, b, c) = map(rec, (in1, in2, in3)),
                ref = :($(rec(out))[1])
                :($qasm_cgc([single_qubit_circuit_gate(($ref), RzGate($a), $N[]),
                   single_qubit_circuit_gate(($ref), RyGate($b), $N[]),
                   single_qubit_circuit_gate(($ref), RzGate($c), $N[]),]))
            end

        Struct_x(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), X, $N[])]))
            end

        Struct_y(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), Y, $N[])]))
            end

        Struct_z(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), Z, $N[])]))
            end

        Struct_h(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), HadamardGate(), $N[])]))
            end

        Struct_t(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), TGate(), $N[])]))
            end

        Struct_tdg(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), TdagGate(), $N[])]))
            end

        Struct_s(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), SGate(), $N[])]))
            end

        Struct_sdg(out=out) =>
            let ref = :($(rec(out)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), SdagGate(), $N[])]))
            end

        Struct_rx(in=in, out=out) =>
            let ref = :($(rec(out))),
                arg = :($(rec(in)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), RxGate($arg), $N[])]))
            end

        Struct_ry(in=in, out=out) =>
            let ref = :($(rec(out))),
                arg = :($(rec(in)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), RyGate($arg), $N[])]))
            end

        Struct_rz(in=in, out=out) =>
            let ref = :($(rec(out))),
                arg = :($(rec(in)))
                :($qasm_cgc([single_qubit_circuit_gate(($ref), RzGate($arg), $N[])]))
            end

        Struct_crz(in=in, out1=out1, out2=out2) =>
            let out = :($(rec(out2))),
                cntrl = :($(rec(out1))),
                arg = :($(rec(in)))

                :($qasm_cgc([controlled_circuit_gate(($out), ($cntrl), RzGate($arg), $N[])]))
            end

        Struct_gate(
            decl = Struct_gatedecl(
                id=Token(str=fid),
                args=args,
                outs=outs
            ),
            goplist=goplist
         ) =>
            let out_ids :: Vector{Symbol} = rec(outs),
                fid = Symbol("custom_gate_"*fid),
                args = rec(args),
                goplist = trans_gates.(goplist, Ref(qasm_cgc), Ref(N))

                if isnothing(goplist)
                    goplist=[]
                end
                if isnothing(args)
                    args=[]
                end
                quote
                    function $fid(($(args...), ), $(out_ids...); qasm_cgc=$qasm_cgc, N=$N)
                        $(goplist...)
                    end
                end
            end

        Struct_ifstmt(l=Token(str=l), r=r, gate_name=Token(str=gate_name), args=nothing, outs=outs) =>
            let l = Symbol(l),
                r = rec(r),
                refs = rec(outs),
                gate_name = Symbol("custom_gate_"*gate_name)

                :($gate_name(nothing, reg_check($l,$r), $(refs...); qasm_cgc=$qasm_cgc, N=$N))
            end

        Struct_ifstmt(l=Token(str=l), r=r, gate_name=Token(str=gate_name), args=exprlist, outs=outs) =>
            let l = Symbol(l),
                r = rec(r),
                refs = rec(outs),
                exprs = Expr(:tuple, rec(exprlist)...),
                gate_name = Symbol("custom_gate_"*gate_name)

                :($gate_name($exprs, reg_check($l,$r), $(refs...); qasm_cgc=$qasm_cgc, N=$N))
            end

        Struct_mainprogram(
            prog = stmts
        ) =>
            let stmts = trans_gates.(stmts, qasm_cgc, N)
                stmts
            end
        _ => nothing
    end
end


function transform_qasm(ctx_tokens)

    cregs = Qaintessent.CRegister[]
    qregs = Qaintessent.QRegister[]
    reg_declr = trans_reg(ctx_tokens, Ref(cregs), Ref(qregs))
    eval.(reg_declr)

    qasm_cgc = CircuitGateChain(qregs, cregs)
    N = size(qasm_cgc)

    gates_declr = trans_gates(ctx_tokens, Ref(qasm_cgc), Ref(N))
    eval.(gates_declr)

    qasm_cgc
end

function qasm2cgc(txt::String)
    qasmlex = lex(txt)
    qasmparse = parse_qasm(qasmlex)
    qasm_cgc = transform_qasm(qasmparse)
end
