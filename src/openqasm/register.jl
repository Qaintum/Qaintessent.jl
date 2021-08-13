abstract type Register end

mutable struct CircuitReference
    reference::Ref{Circuit}
    function CircuitReference()
        new()
    end
end

"""
    QRegister{N}

Quantum Register of `N` qubits.

"""
struct QRegister <: Register
    # TODO: Create function to read QRegister using ind
    "size of register"
    N::Int64
    "vector of quantum wire numbers in this quantum register"
    ind::Vector{Int}
    "ref to circuit"
    circuit::CircuitReference

    # function QRegister(N::Int64, ind::Vector{Int})
    #     new(N, ind)
    # end
end

function (qreg::QRegister)(indices::Int...) 
    qreg, qreg[[indices...]]
end

function (qreg::QRegister)(indices::Array{Int,1}) 
    qreg, qreg[indices]
end

# TODO: Reimplement CRegister smartly
# """
#     CRegister{N}

# Classical register of `N` bits.
# """
# struct CRegister <: Register
#     # TODO: Create function to read QRegister using ind
#     "size of register"
#     N::Int
#     "vector of quantum wire numbers in this quantum register"
#     ind::Vector{Int}
#     "bit array of bits in this classical register"
#     n::BitArray{1}
#     "ref to circuit"
#     circuit::CircuitReference

#     # function CRegister(N::Int64, ind::Vector{Int}, n::BitArray{1})
#     #     new(N, ind, n)
#     # end
# end

"""
    qreg(n::Int)

constructs quantum register of size `n`
"""
function qreg(n::Int)
    QRegister(n, zeros(Int, n), CircuitReference())
end

# """
#     creg(n::Int)

# constructs classical register of size `n`
# """
# function creg(n::Int)
#     CRegister(n, zeros(Int, n), BitArray(undef, (n)), CircuitReference())
# end



# make Registers indexable
function Base.getindex(r::Register, i::Int)
    !any(r.ind .== 0) || error("Register object has yet to be used in a Circuit object")
    minimum(i) >= 1 || throw(BoundsError(r.ind, minimum(i)))
    maximum(i) <= length(r.ind) || throw(BoundsError(r.ind, maximum(i)))

    return r.ind[i]
end

function Base.getindex(r::Register, i::Array{Int,1})
    !any(r.ind .== 0) || error("Register object has yet to be used in a Circuit object")
    minimum(i) >= 1 || throw(BoundsError(r.ind, minimum(i)))
    maximum(i) <= length(r.ind) || throw(BoundsError(r.ind, maximum(i)))

    return r.ind[i]
end

function Base.getindex(r::Register, ::Colon)
    !any(r.ind .== 0) || error("Register object has yet to be used in a Circuit object")

    return r.ind[:]
end

function Base.getindex(r::Register, ur::UnitRange{Int64})
    !any(r.ind .== 0) || error("Register object has yet to be used in a Circuit object")

    return r.ind[ur]
end


# """
#     reg_check(c::Qaintessent.CRegister, val::Int)

# create expression to verify that classical register is equal to integer value `val`
# """
# function reg_check(c::Qaintessent.CRegister, val::Int)
#     r = Ref(c.n)
#     return :(Qaintessent.binary_to_int($r[]) == $val)
# end

# """
#     set_creg!(reg::CRegister, index::Int, value::Bool) where {N}

# Set the classical register `reg` bit `index` to value `value`
# """
# function set_creg!(reg::CRegister, index::Int, value::Bool)
#     index = length(reg.n) - index + 1
#     reg.n[index] = value
# end

# """
#     set_creg!(reg::CRegister, value::Int) where {N}

# Set the classical register `reg` to value `value`
# """
# function set_creg!(reg::CRegister, value::Int)
#     l = length(reg.n)
#     log2(value) <= l || error("Unable to store integer value " * string(value) * " in BitArray of size " * string(l))
#     reg.n .= binary_digits(l, value)
# end

# num_wires(creg::CRegister) = creg.N
num_wires(qreg::QRegister) = qreg.N

"""
    Circuit(regs::Register...)

helper function, creates Circuit from multiple quantum or classical registers
"""
function Circuit(regs::Register...)
    N = 0
    j = 1
    for i in 1:length(regs)
        !isdefined(regs[i].circuit, :reference) || error("Input Register $i is already used in another circuit.")
        N += regs[i].N
        regs[i].ind .= j:regs[i].N+j-1
        j = regs[i].N+j
    end
    c = Circuit{N}()
    for i in 1:length(regs)
        regs[i].circuit.reference = c
    end
    return c
end

function add!(g::AbstractGate, q::QRegister)
    isdefined(q.circuit, :reference) || error("Register object has yet to be used in a Circuit object")
    if num_wires(g) == num_wires(q)
        append!(q.circuit.reference[], circuit_gate(Tuple(q[:]), g))
        return
    end 
    num_wires(g) == 1 || error("AbstractGate $g is applied to $(num_wires(g)) wires, quantum register of size $(num_wires(q)) provided. Unable to determine gates to be applied.")

    for i in 1:num_wires(q)
        append!(q.circuit.reference[], circuit_gate(q[i], g))
    end
end

function add!(g::AbstractGate, target_wires::Tuple{QRegister, Array{Int}}...)
    circuit_length = num_wires(g)
    circuit_references = getfield.(target_wires, 1)

    wires = reduce(hcat, getfield.(target_wires, 2))
    length(wires) == circuit_length || error("AbstractGate $g is applied to $circuit_length wires. However, $(length(wires)) wires provided")
    for ref in circuit_references[2:end]
        isdefined(ref.circuit, :reference) || error("Register object has yet to be used in a Circuit object")
        ref.circuit.reference[] === circuit_references[1].circuit.reference[] || error("Targetted Registers belong to different Circuit objects")
    end
    append!(circuit_references[1].circuit.reference[], CircuitGate(Tuple(wires), g))
    return
end

function add_control!(g::AbstractGate, target_q::QRegister, control_q::QRegister)
    num_wires(g) == 1 || error("Helper function for controlled gates across entire registers only supports 1 qubit gates")
    isdefined(target_q.circuit, :reference) || error("Target quantum register is not included in any circuit")
    isdefined(control_q.circuit, :reference) || error("Control quantum register is not included in any circuit")
    target_q.circuit.reference[] === control_q.circuit.reference[] || error("Target quantum register and control quantum register are used in different circuits")

    num_wires(target_q) == num_wires(control_q) || error("Control register and Target register have differing number of qubits, unable to determine gates to be applied")
    for i in 1:num_wires(control_q)
        append!(target_q.circuit.reference[], circuit_gate(target_q.ind[i], g, control_q.ind[i]))
    end
end

function add_control!(g::AbstractGate, target_wires::Tuple{QRegister, Array{Int}}, control_wires::Tuple{QRegister, Array{Int}})
    isdefined(target_wires[1].circuit, :reference) || error("Target quantum register is not included in any circuit")
    isdefined(control_wires[1].circuit, :reference) || error("Control quantum register is not included in any circuit")
    target_wires[1].circuit.reference[] === control_wires[1].circuit.reference[] || error("Target quantum register and control quantum register are used in different circuits")

    circuit_length = num_wires(g)
    num_wires(g) == length(target_wires[2]) || error("Gate affecting $(num_wires(g)) qubits applied to $(length(target_wires[2])) qubits")
    append!(target_wires[1].circuit.reference[], circuit_gate((target_wires[2]...,), g, control_wires[2]...))
end