"""
    Moment

Represents an intermediate state within a given circuit.
"""
mutable struct Moment
    gates::Vector{CircuitGate}

    @doc """
        Moment(g::AbstractCircuitGate) where 

    Create a `Moment` object consisting of a single `CircuitGate` object.
    """
    function Moment(g::AbstractCircuitGate)
        new([g])
    end

    @doc """
        Moment(g::AbstractVector{<:AbstractCircuitGate}) where 

    Create a `Moment` object consisting of multiple `CircuitGate` objects.
    """
    function Moment(cgs::Vector{<:AbstractCircuitGate})
        iwire = Int[]
        for cg in cgs
            length(intersect(iwire, cg.iwire)) == 0 || error("Only gates on different wires are allowed in a Moment")
            append!(iwire, collect(cg.iwire))
        end
        new(cgs)
    end
end

"""
    Base.adjoint(m::Moment) 

Construct a `Moment` object that is the adjoint of `m`.
"""
function Base.adjoint(m::Moment) 
    return Moment(Base.adjoint.(reverse(m.gates)))
end

"""
    sparse_matrix(m::Moment) 

returns matrix representation of a `Moment{M}` object that can applied to a state vector of `N` qubits.
"""
function sparse_matrix(m::Moment, N::Int=0)
    length(m) != 0 || error("Moment with no gates cannot be converted to matrix")
    if N == 0
        N = req_wires(m)
    end
    mat = sparse_matrix(m[1], N)
    for i in 2:length(m)
        mat = sparse_matrix(m[i], N) * mat
    end
    mat
end

function sparse_matrix(m::Vector{Moment}, N::Int=0)::SparseMatrixCSC{Complex{FloatQ},Int}
    length(m) != 0 || error("Vector of length 0 cannot be converted to matrix")
    if N == 0
        N = maximum(req_wires.(m))
    end
    mat = sparse_matrix(m[1], N)
    for i in 2:length(m)
        mat = sparse_matrix(m[i], N) * mat
    end
    mat
end

# make Moment iterable and indexable
function Base.getindex(m::Moment, i::Integer) 
    1 <= i <= length(m.gates) || throw(BoundsError(m, i))
    return m.gates[i]
end

function Base.getindex(m::Moment, ::Colon) 
    return m.gates[:]
end

function Base.getindex(m::Moment, ur::UnitRange{Int64}) 
    return m.gates[ur]
end

function Base.iterate(m::Moment, state=1) 
    return state > length(m.gates) ? nothing : (m[state], state + 1)
end

# implement methods required for iteration
function Base.firstindex(::Moment) 
    return 1
end

function Base.lastindex(m::Moment) 
    return length(m.gates)
end

function Base.length(m::Moment) 
    return length(m.gates)
end

function Base.isapprox(m1::Moment, m2::Moment)
    length(m1) == length(m2) || return false
    for (g1, g2) in zip(m1, m2)
        if !(g1 ≈ g2)
            return false
        end
    end
    return true
end

function Base.reverse!(m::Moment) 
    m.gates = reverse(m.gates)
    return m
end

function Base.reverse(m::Moment) 
    Moment(reverse(m.gates))
end

@memoize function req_wires(m::Moment)
    N = 0
    for cg in m
        N = maximum((N, req_wires(cg)))
    end
    return N
end

function Base.pop!(m::Moment)
    pop!(m.gates)
end

function Base.isempty(m::Moment)
    isempty(m.gates)
end