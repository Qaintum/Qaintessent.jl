using RBNF
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
@as_record Struct_idlist
@as_record Struct_mixedlist
@as_record Struct_argument
@as_record Struct_explist
@as_record Struct_pi
@as_record Struct_fnexp
@as_record Struct_neg
@as_record Struct_bin


struct Session
    n       :: Ref{Int}
    gates   :: Bool
end


# init_sess() = Session(Ref(1), false)
# trans(qasm) = trans(qasm, init_sess())

function trans(ctx_tokens)

    function app(op, args...)
        args = map(rec, args)
        op = Symbol(op)
        :($op($(args...)))
    end
    wires = Dict{Symbol, Int}()

    rec(qasm) = trans(qasm)
    cregs = CRegister[]
    qregs = QRegister[]
    @match token begin
        Struct_decl(
            regtype = Token(str=regtype),
            id = Token(str=id),
            int = Token(str = n)
        ) =>
            let id = Symbol(id),
                n = parse(Int, n)

                if regtype == "qreg"
                    :($id = qreg(n);push!(qregs, $id))
                elseif regtype == "creg"
                    :($id = creg(n);push!(cregs, $id))
                end
            end

        Struct_mainprogram(prog=prog)
            => let prog = map(rec, prog)
                Expr(:block, prog...)
                end
        _ => println("Not included")
    end

    cgc = CircuitGateChain(qregs, cregs)
    println(cgc)

    # @match token begin
    #     Struct_pi(_) => Base.pi
    #     Token{:id}(str=str) => Symbol(str)
    #     Token{:real}(str=str) => parse(Float64, str)
    #     Token{:nninteger}(str=str) => parse(Int64, str)
    #     Struct_neg(value=value) => :(-$(rec(value)))
    #     Struct_bin(l=l, op=Token(str=op), r=r) => app(op, l, r)
    #     Struct_fnexp(fn = Token(str=fn), arg=arg) =>
    #         let fn = @match fn begin
    #             "sin" => sin
    #             "cos" => cos
    #             "tan" => tan
    #             "exp" => exp
    #             "ln"  => log
    #             "sqrt"=> sqrt
    #             _     => error("not impl yet")
    #             end
    #             app(fn, arg)
    #         end
    #
    #     Struct_idlist(hd=Token(str=hd), tl=nothing) => [Symbol(hd)]
    #     Struct_idlist(hd=Token(str=hd), tl=tl) => [Symbol(hd), rec(tl)...]
    #
    #     Struct_explist(hd=hd, tl=nothing) => [rec(hd)]
    #     Struct_explist(hd=hd, tl=tl) => [rec(hd), rec(tl)...]
    #
    #     Struct_mixedlist(hd=hd, tl=nothing) => [rec(hd)]
    #     Struct_mixedlist(hd=hd, tl=tl) => [rec(hd), rec(tl)...]
    #
    #     Struct_argument(id=Token(str=id), arg=nothing) => Symbol(id)
    #     Struct_argument(id=Token(str=id), arg=Token(str=int)) =>
    #         let ind = parse(Int, int) + 1 # due to julia 1-based index
    #             :($(Symbol(id))[$ind])
    #         end
    #     Struct_mainprogram(prog=prog)
    #         => let prog = map(rec, prog)
    #             Expr(:block, prog...)
    #         end
    #     _ => println("Not included")
    #     end
    # end
end
