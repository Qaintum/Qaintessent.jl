abstract type Register end

struct QRegister <: Register
    # register size
    n::Int
    ind::Vector{Int}
end

"""
    CRegister{N}

Classical register of N bits.
"""
struct CRegister <: Register
    # register size
    n::BitArray{1}
    # bit array representing register
    ind::Vector{Int}
end

function qreg(n::Int)
    QRegister(n, zeros(n))
end

function creg(n::Int)
    CRegister(BitArray(undef, (n)), zeros(n))
end

# make Registers indexable
function Base.getindex(r::Register, i::Integer) where {N}
    !any(r.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    1 <= i <= length(r.ind) || throw(BoundsError(ind, i))
    return r.ind[i]
end

"""
    reg_check(c::Qaintessent.CRegister, val::Int)

create expression to verify that classical register is equal to integer value `val`
"""
function reg_check(c::Qaintessent.CRegister, val::Int)
    r = Ref(c.n)
    return :(bit2int($r[])==$val)
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
