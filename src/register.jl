abstract type Register end

struct QRegister <: Register
    n::BitArray{1}
    ind::Vector{Int}
end

struct CRegister <: Register
    n::BitArray{1}
    ind::Vector{Int}
end

function qreg(n::Int)
    QRegister(BitArray(undef, (n)), zeros(n))
end

function creg(n::Int)
    CRegister(BitArray(undef, (n)), zeros(n))
end

# make CircuitGateChain iterable and indexable
function Base.getindex(r::Register, i::Integer) where {N}
    !any(r.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    1 <= i <= length(r.ind) || throw(BoundsError(ind, i))
    return r.ind[i]
end
