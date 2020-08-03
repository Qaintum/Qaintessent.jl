using RBNF
using MLStyle
using PrettyPrint
struct QASMLang end

second((a, b)) = b
second(vec::V) where V <: AbstractArray = vec[2]

struct Struct_bin
    l
    op :: RBNF.Token
    r
end

RBNF.@parser QASMLang begin
    # define ignorances
    ignore{includes, space, comments}

    @grammar
    # define grammars
    mainprogram := ["OPENQASM", ver=real, ';', prog=program]
    program     = statement{*}
    statement   = (decl | gate | opaque | ifstmt | barrier | qop)
    # stmts
    ifstmt      := ["if", '(', l=id, "==", r=nninteger, ')', gate_name=id, ['(', [args=explist].?, ')'].?, outs=mixedlist, ';']
    opaque      := ["opaque", id=id, ['(', [arglist1=idlist].?, ')'].? , arglist2=idlist, ';']
    barrier     := ["barrier", value=mixedlist, ';']
    decl        := [regtype=("qreg" | "creg"), id=id, '[', int=nninteger, ']', ';']

    # gate
    gate        := [decl=gatedecl, [goplist=goplist].?, '}']
    gatedecl    := ["gate", id=id, ['(', [args=idlist].?, ')'].?, (outs=idlist), '{']

    goplist     = (barrier_ids | uop){*}
    barrier_ids := ["barrier", ids=idlist, ';'] # not impl
    # qop
    qop         = (measure | reset | uop)
    reset       := ["reset", arg=argument, ';'] # not impl
    measure     := ["measure", arg1=argument, "->", arg2=argument, ';'] # not impl

    uop         = (u | cx | h | x | y | z | s | sdg | t | tdg | rx | ry | rz | crz | iduop)
    u          := ['U', '(', in1=exp, ',', in2=exp, ',', in3=exp, ')', out=argument, ';']
    cx         := ["CX", out1=argument, ',', out2=argument, ';']
    ch         := ["ch", out1=argument, ',', out2=argument, ';']
    h          := ['h',  out=argument, ';']
    x          := ['x',  out=argument, ';']
    y          := ['y',  out=argument, ';']
    z          := ['z',  out=argument, ';']
    s          := ['s',  out=argument, ';']
    sdg        := ["sdg",  out=argument, ';']
    t          := ['t',  out=argument, ';']
    tdg        := ["tdg",  out=argument, ';']
    rx         := ["rx",  '(', in=exp, ')', out=argument, ';']
    ry         := ["ry",  '(', in=exp, ')', out=argument, ';']
    rz         := ["rz",  '(', in=exp, ')', out=argument, ';']
    crz        := ["crz",  '(', in=exp, ')', out1=argument, out2=argument, ';']

    iduop      := [gate_name=id, ['(', [args=explist].?, ')'].?, outs=mixedlist, ';']

    idlist     := [hd=id, [',', tl=idlist].?]

    mixedlist  := [hd=argument, [',', tl=mixedlist].?]

    argument   := [id=id, ['[', (arg=nninteger), ']'].?]

    explist    := [hd=exp, [',', tl=explist].?]
    pi         := "pi"
    atom       =  (real | nninteger | fnexp | pi | id) | (['(', exp, ')'] % second) | neg
    fnexp      := [fn=fn, '(', arg=exp, ')']
    neg        := ['-', value=exp]

    exp        = [l=mul,  [op=('+' |'-'), r=exp].?] => _.op === nothing ? _.l : Struct_bin(_.l, _.op, _.r)
    mul        = [l=atom, [op=('*' | '/'), r=mul].?] => _.op === nothing ? _.l : Struct_bin(_.l, _.op, _.r)
    fn         = ("sin" | "cos" | "tan" | "exp" | "ln" | "sqrt")

    # define tokens
    @token
    includes  := r"\Ginclude .*;"
    id        := r"\G[a-z]{1}[A-Za-z0-9_]*"
    real      := r"\G([0-9]+\.[0-9]*|[0-9]*\.[0.9]+)([eE][-+]?[0-9]+)?"
    nninteger := r"\G([1-9]+[0-9]*|0)"
    space     := r"\G\s+"
    comments  := r"\G//.*"
end


Token = RBNF.Token


function lex(src :: String)
    RBNF.runlexer(QASMLang, src)
end

function parse_qasm(tokens :: Vector{Token{A} where A})
    ast, ctx = RBNF.runparser(mainprogram, tokens)
    ast
end
