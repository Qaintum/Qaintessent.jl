using Test
using TestSetExtensions
using RBNF
using Qaintessent
using JLD2


# Adding random comments and spacings in src1 to test robustness of lexer
f = open("qpe.txt")
s = read(f, String)

# println(Qaintessent.parse_qasm(Qaintessent.lex(s)))

cgc = qasm2cgc(s)
JLD2.save_object("cgc.jld2", cgc)
# println(cgc)
