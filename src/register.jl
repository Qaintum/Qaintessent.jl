abstract type Register end

"""
    QRegister{N}

Quantum Register of `N` qubits.

"""
struct QRegister <: Register
    # TODO: Create function to read QRegister using ind
    "value stored by register"
    n::Int
    "vector of quantum wire numbers in this quantum register"
    ind::Vector{Int}
end

"""
    CRegister{N}

Classical register of `N` bits.
"""
struct CRegister <: Register
    "bit array of bits in this classical register"
    n::BitArray{1}
    "vector of quantum wire numbers in this quantum register"
    ind::Vector{Int}
end

"""
    qreg(n::Int)

constructs quantum register of size `n`
"""
function qreg(n::Int)
    QRegister(n, zeros(n))
end

"""
    creg(n::Int)

constructs classical register of size `n`
"""
function creg(n::Int)
    CRegister(BitArray(undef, (n)), zeros(n))
end

# make Registers indexable
function Base.getindex(r::Register, i::Integer) where {N}
    !any(r.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    1 <= i <= length(r.ind) || throw(BoundsError(ind, i))
    return r.ind[i]
end


function int2bit(x::Int64; pad=nothing)
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
function set_creg!(reg::CRegister, index::Int, value::Bool) where {N}
    index = length(reg.n) - index + 1
    reg.n[index] = value
end

"""
    set_creg!(reg::CRegister, value::Int) where {N}

Set the classical register `reg` to value `value`
"""
function set_creg!(reg::CRegister, value::Int) where {N}
    l = length(reg.n)
    log2(value) <= l || error("Unable to store integer value " * string(value) * " in BitArray of size " * string(l))
    reg.n .= int2bit(value; pad=l)
end
