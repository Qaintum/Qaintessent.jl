using Qaintessent
using Random
using LinearAlgebra
using JLD2
using BenchmarkTools


# Adding random comments and spacings in src1 to test robustness of lexer
# f = open("qpe.txt")
# s = read(f, String)

# println(Qaintessent.parse_qasm(Qaintessent.lex(s)))

# cgc = qasm2cgc(s)
# JLD2.save_object("cgc.jld2", cgc)
cgc = JLD2.load_object("cgc.jld2")

ψ = rand(ComplexF64, 2^16)
ψ = ψ ./ norm(ψ)
ψ = Statevector(ψ)
apply!(ψ, cgc.moments)