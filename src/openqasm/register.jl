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

"""
    CRegister{N}

Classical register of `N` bits.
"""
struct CRegister <: Register
    # TODO: Create function to read QRegister using ind
    "size of register"
    N::Int
    "vector of quantum wire numbers in this quantum register"
    ind::Vector{Int}
    "bit array of bits in this classical register"
    n::BitArray{1}
    "ref to circuit"
    circuit::CircuitReference

    # function CRegister(N::Int64, ind::Vector{Int}, n::BitArray{1})
    #     new(N, ind, n)
    # end
end

"""
    qreg(n::Int)

constructs quantum register of size `n`
"""
function qreg(n::Int)
    QRegister(n, zeros(Int, n), CircuitReference())
end

"""
    creg(n::Int)

constructs classical register of size `n`
"""
function creg(n::Int)
    CRegister(n, zeros(Int, n), BitArray(undef, (n)), CircuitReference())
end



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

function int2bit(x::Int; pad=nothing)
    if isnothing(pad)
        pad = floor(Int64, log2(3) + 1)
    end
    BitArray(digits(x, base=2, pad=pad) |> reverse)
end

function bit2int(b::BitArray)
    b = deepcopy(b)
    num = 0
    count = 0
    while !isempty(b)
        num += pop!(b) * 2^count
        count += 1
    end
    num
end


"""
    reg_check(c::Qaintessent.CRegister, val::Int)

create expression to verify that classical register is equal to integer value `val`
"""
function reg_check(c::Qaintessent.CRegister, val::Int)
    r = Ref(c.n)
    return :(Qaintessent.bit2int($r[]) == $val)
end

"""
    set_creg!(reg::CRegister, index::Int, value::Bool) where {N}

Set the classical register `reg` bit `index` to value `value`
"""
function set_creg!(reg::CRegister, index::Int, value::Bool)
    index = length(reg.n) - index + 1
    reg.n[index] = value
end

"""
    set_creg!(reg::CRegister, value::Int) where {N}

Set the classical register `reg` to value `value`
"""
function set_creg!(reg::CRegister, value::Int)
    l = length(reg.n)
    log2(value) <= l || error("Unable to store integer value " * string(value) * " in BitArray of size " * string(l))
    reg.n .= int2bit(value; pad=l)
end

Base.size(creg::CRegister) = creg.N
Base.size(qreg::QRegister) = qreg.N

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
    c = Circuit{N}(Moment[])
    for i in 1:length(regs)
        regs[i].circuit.reference = c
    end
    return c
end

function add!(g::AbstractGate, q::QRegister)
    isdefined(q.circuit, :reference) || error("Register object has yet to be used in a Circuit object")
    if num_wires(g) == size(q)
        append!(q.circuit.reference[], circuit_gate(q[:], g))
    end 
    num_wires(g) == 1 || error("AbstractGate $g is applied to $(num_wires(g)) wires, quantum register of size $(size(q)) provided. Unable to determine gates to be applied.")

    for i in 1:size(q)
        append!(q.circuit.reference[], circuit_gate(q[i], g))
    end
end

function add!(g::AbstractGate, target_wires::Tuple{QRegister, Array{Int}}...)
    circuit_length = num_wires(g)
    circuit_references = getfield.(target_wires, 1)
    wires = getfield.(target_wires, 2)
    for ref in circuit_references[2:end]
        isdefined(ref.circuit, :reference) || error("Register object has yet to be used in a Circuit object")
        ref.circuit.reference[] === circuit_references[1].circuit.reference[] || error("Targetted Registers belong to different Circuit objects")
    end
    for i in 1:circuit_length
        append!(circuit_references[1].circuit.reference[], circuit_gate(getindex.(wires, i), g))
    end
end

function add_control!(g::AbstractGate, target_q::QRegister, control_q::QRegister)
    isdefined(target_q.circuit, :reference) || error("Target quantum register is not included in any circuit")
    isdefined(control_q.circuit, :reference) || error("Control quantum register is not included in any circuit")
    target_q.circuit.reference[] === control_q.circuit.reference[] || error("Target quantum register and control quantum register are used in different circuits")

    circuit_length = num_wires(g)
    num_wires(g) == size(target_q) || error("Gate affecting $(num_wires(g)) qubits applied to quantum register of $(size(target_q)) qubits")
    append!(target_q.circuit.reference[], circuit_gate((target_q.ind...), g, control_q.ind...))
end

function add_control!(g::AbstractGate, target_wires::Tuple{QRegister, Array{Int}}, control_wires::Tuple{QRegister, Array{Int}})
    isdefined(target_wires[1].circuit, :reference) || error("Target quantum register is not included in any circuit")
    isdefined(control_wires[1].circuit, :reference) || error("Control quantum register is not included in any circuit")
    target_wires[1].circuit.reference[] === control_wires[1].circuit.reference[] || error("Target quantum register and control quantum register are used in different circuits")

    circuit_length = num_wires(g)
    num_wires(g) == length(target_wires[2]) || error("Gate affecting $(num_wires(g)) qubits applied to $(length(target_wires[2])) qubits")
    append!(target_wires[1].circuit.reference[], circuit_gate((target_wires[2]...,), g, control_wires[2]...))
end