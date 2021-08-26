using RBNF
using RBNF: Token
using MLStyle
using PrettyPrint

@as_record Struct_mainprogram
@as_record Struct_ifstmt
@as_record Struct_gate
@as_record Struct_gatedecl
@as_record Struct_decl
@as_record Struct_barrier
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
@as_record Struct_crx
@as_record Struct_cry
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

function trans_reg(ctx_tokens, qregs)
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
                    return nothing
                end
            end

        Struct_mainprogram(
            prog = stmts
        ) =>
            let stmts = trans_reg.(stmts, qregs)
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
                :(append!($qasm_cgc, [circuit_gate(($ref2), X, ($ref1))]))
            end

        Struct_ch(out1=out1, out2=out2) =>
            let ref1 = rec(out1),
                ref2 = rec(out2)
                :(append!($qasm_cgc, [circuit_gate(($ref2), HadamardGate(), ($ref1))]))
            end

        Struct_u(in1=in1, in2=in2, in3=in3, out=out) =>
            let (a, b, c) = map(rec, (in1, in2, in3)),
                ref = :($(rec(out))[1])
                :(append!($qasm_cgc, [circuit_gate(($ref), RzGate($a)),
                   circuit_gate(($ref), RyGate($b)),
                   circuit_gate(($ref), RzGate($c)),]))
            end

        Struct_x(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), X)]))
            end

        Struct_y(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), Y)]))
            end

        Struct_z(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), Z)]))
            end

        Struct_h(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), HadamardGate())]))
            end

        Struct_t(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), TGate())]))
            end

        Struct_tdg(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), TdagGate())]))
            end

        Struct_s(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), SGate())]))
            end

        Struct_sdg(out=out) =>
            let ref = :($(rec(out)))
                :(append!($qasm_cgc, [circuit_gate(($ref), SdagGate())]))
            end

        Struct_rx(in=in, out=out) =>
            let ref = :($(rec(out))),
                arg = :($(rec(in)))
                :(append!($qasm_cgc, [circuit_gate(($ref), RxGate($arg))]))
            end

        Struct_ry(in=in, out=out) =>
            let ref = :($(rec(out))),
                arg = :($(rec(in)))
                :(append!($qasm_cgc, [circuit_gate(($ref), RyGate($arg))]))
            end

        Struct_rz(in=in, out=out) =>
            let ref = :($(rec(out))),
                arg = :($(rec(in)))
                :(append!($qasm_cgc, [circuit_gate(($ref), RzGate($arg))]))
            end
            
        Struct_crx(in=in, out1=out1, out2=out2) =>
            let out = :($(rec(out2))),
                cntrl = :($(rec(out1))),
                arg = :($(rec(in)))

                :(append!($qasm_cgc, [circuit_gate(($out), RxGate($arg), ($cntrl))]))
            end

        Struct_cry(in=in, out1=out1, out2=out2) =>
            let out = :($(rec(out2))),
                cntrl = :($(rec(out1))),
                arg = :($(rec(in)))

                :(append!($qasm_cgc, [circuit_gate(($out), RyGate($arg), ($cntrl))]))
            end

        Struct_crz(in=in, out1=out1, out2=out2) =>
            let out = :($(rec(out2))),
                cntrl = :($(rec(out1))),
                arg = :($(rec(in)))

                :(append!($qasm_cgc, [circuit_gate(($out), RzGate($arg), ($cntrl))]))
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

        # Struct_ifstmt(l=Token(str=l), r=r, gate_name=Token(str=gate_name), args=nothing, outs=outs) =>
        #     let l = Symbol(l),
        #         r = rec(r),
        #         refs = rec(outs),
        #         gate_name = Symbol("custom_gate_"*gate_name)

        #         :($gate_name(nothing, reg_check($l,$r), $(refs...); qasm_cgc=$qasm_cgc, N=$N))
        #     end

        # Struct_ifstmt(l=Token(str=l), r=r, gate_name=Token(str=gate_name), args=exprlist, outs=outs) =>
        #     let l = Symbol(l),
        #         r = rec(r),
        #         refs = rec(outs),
        #         exprs = Expr(:tuple, rec(exprlist)...),
        #         gate_name = Symbol("custom_gate_"*gate_name)

        #         :($gate_name($exprs, reg_check($l,$r), $(refs...); qasm_cgc=$qasm_cgc, N=$N))
        #     end

        Struct_mainprogram(
            prog = stmts
        ) =>
            let stmts = trans_gates.(stmts, qasm_cgc, N)
                stmts
            end
        _ => nothing
    end
end

function trans_measure(ctx_tokens, meas, N)
    function app(op, args...)
        args = map(rec, args)
        op = Symbol(op)
        :($op($(args...)))
    end

    @match ctx_tokens begin
        Struct_measure(arg1=arg1, arg2=arg2) =>
            let ref = :($(rec(arg1)))
                :(append!($meas, [mop(Z, $ref)]))  
            end

        Struct_barrier(value=value) =>
            let ref = :($(rec(value)))
                :(@warn "Barrier operation is not yet supported")
            end

        Struct_reset(arg=arg) =>
            let ref = :($(rec(arg)))
                :(@warn "Reset operation is not yet supported")
            end

        Struct_mainprogram(
            prog = stmts
        ) =>
            let stmts = trans_measure.(stmts, meas, N)
                stmts
            end
        _ => nothing
    end
end



function transform_qasm(ctx_tokens)

    qregs = Qaintessent.QRegister[]
    reg_declr = trans_reg(ctx_tokens, Ref(qregs))
    eval.(reg_declr)

    qasm_cgc = Circuit(qregs...)
    meas = MeasurementOperator[]
    N = num_wires(qasm_cgc)

    gates_declr = trans_gates(ctx_tokens, Ref(qasm_cgc), Ref(N))

    eval.(gates_declr)

    measure_declr = trans_measure(ctx_tokens, Ref(meas), Ref(N))
    eval.(measure_declr)
    add_measurement!(qasm_cgc, meas)
    
    qasm_cgc
end

"""
    qasm2cgc(txt::String)

converts OpenQASM 2.0 text to Circuit{N} object
"""
function qasm2cgc(txt::String)
    qasmlex = lex(txt)
    qasmparse = parse_qasm(qasmlex)
    transform_qasm(qasmparse)
end


function Base.:(==)(dcl1::Struct_gate, dcl2::Struct_gate)
    dcl1.decl == dcl2.decl && dcl1.goplist == dcl2.goplist
end
