using Memoize
using LinearAlgebra

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
    function Moment(gs::Vector{<:AbstractCircuitGate})
        iwire = Integer[]
        Nmin = 0
        for gate in gs
            length(intersect(iwire, gate.iwire)) == 0 || error("Only gates on different wires are allowed in a Moment")
            Nmin = maximum((Nmin, num_wires(gate.gate)))
            append!(iwire, collect(gate.iwire))
        end
        new(gs)
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
    matrix(m::Moment) 

returns matrix representation of a `Moment{M}` object that can applied to a state vector of `N` qubits.
"""
function matrix(m::Moment, N::Int=0)
    if N == 0
        N = size(m)
    end
    mat = matrix(m[1], N)
    for i in 2:length(m)
        mat = matrix(m[i], N) * mat
    end
    mat
end

function matrix(m::Vector{Moment}, N::Int=0)::SparseMatrixCSC{Complex{Float64},Int}
    if N == 0
        N = maximum(size.(m))
    end
    mat = matrix(m[1], N)
    for i in 2:length(m)
        mat = matrix(m[i], N) * mat
    end
    mat
end

# make Moment iterable and indexable
function Base.getindex(m::Moment, i::Integer) 
    1 <= i <= length(m.gates) || throw(BoundsError(m, i))
    return m.gates[i]
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

@memoize function Base.size(m::Moment)
    Nmin = 0
    for circuit_gate in m
        Nmin = maximum((Nmin, maximum(num_wires(circuit_gate))))
    end
    Nmin
end


"""
    rdm(N, iwire, ψ, χ)

Compute the reduced density matrix ``tr_B[|ψ⟩⟨χ|]``, where the trace runs over
the subsystem complementary to the qubits specified by `iwire`.
"""
function rdm(N::Integer, iwire::NTuple{M,<:Integer}, ψ::AbstractVector, χ::AbstractVector) where {M}
    M ≥ 1 || error("Need at least one wire to act on.")
    M ≤ N || error("Number of gate wires cannot be larger than total number of wires.")
    length(unique(iwire)) == M || error("Wire indices must be unique.")
    minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")
    maximum(iwire) ≤ N || error("Wire index cannot be larger than total number of wires.")

    # convert to array
    iwire = collect(iwire)
    # complementary wires
    iwcompl = setdiff(1:N, iwire)
    @assert length(iwire) + length(iwcompl) == N

    # TODO: support general "qudits"
    d = 2

    ρ = zeros(eltype(ψ), d^M, d^M)

    # Note: following the ordering convention of `kron` here, i.e.,
    # last qubit corresponds to fastest varying index
    strides = [d^(j - 1) for j in 1:N]
    wstrides = strides[iwire]
    cstrides = strides[iwcompl]

    # TODO: optimize memory access pattern
    for kw in (N - M > 0 ? cartesian_tuples(d, N - M) : [Int[]])
        koffset = dot(collect(kw), strides[iwcompl])
        for (i, iw) in enumerate(cartesian_tuples(d, M))
            for (j, jw) in enumerate(cartesian_tuples(d, M))
                rowind = koffset + dot(collect(iw), strides[iwire]) + 1
                colind = koffset + dot(collect(jw), strides[iwire]) + 1
                ρ[i, j] += ψ[rowind] * conj(χ[colind])
            end
    end
    end

    return ρ
end


"""
    MeasurementOperator{M,G} <: AbstractCircuitGate

General complex operator. `M` is the number of wires affected by the CircuitGate, `G` is the basic gate used to construct the Operator.
"""
struct MeasurementOperator{M,G}
    operator::G
    iwire::NTuple{M,Int}
    function MeasurementOperator(cg::G, iwire::NTuple{M,Integer}) where {M,G <:AbstractGate}
        M == num_wires(cg) || error("CircuitGate affecting $(num_wires(cg)) wires given, `iwire` of length $(length(iwire)) provided")
        cg == adjoint(cg) || error("Measurement operator must be Hermitian.")
        new{M,G}(cg, iwire)
    end

    function MeasurementOperator(m::G, iwire::NTuple{M,Integer}) where {M,G <:AbstractMatrix}
        d = 2
        size(m) == (d^M, d^M) || error("Measurement operator must be a $(d^M) × $(d^M) matrix.")
        m ≈ Base.adjoint(m) || error("Measurement operator must be Hermitian.")
        new{M,SparseMatrixCSC{Complex{Float64},Int}}(sparse(m), iwire)
    end

    function MeasurementOperator(m, iwire::Integer...)
        MeasurementOperator(m, iwire)
    end
end

mop(operator, iwire) = MeasurementOperator(operator, iwire)

@memoize function Base.size(m::MeasurementOperator)
    length(m.iwire)
end

@memoize function check_commute(ops::Vector{<:MeasurementOperator})
    for (i, operator) in enumerate(ops)
        for operator2 in ops[i+1:end]
            if !iscommuting(operator, operator2)
                return false
            end
        end
    end
    true
end

"""
    matrix(m::MeasurementOperator{M,G}, N::Integer=0) where {M,G}

returns matrix representation of a `MeasurementOperator{M,G}` object that can applied to a state vector of `N` qubits.
"""
function matrix(m::MeasurementOperator{M,G}, N::Integer=0) where {M,G}
    # convert to array
    iwire = Int[i for i in m.iwire]

    if N == 0
        N = maximum([iwire' M])
    else
        N >= maximum(iwire) || error("CircuitGate applied to iwires, $iwire. Input circuit size `N` is $N")
    end

    gmat = matrix(m.operator)

    _matrix(gmat, iwire, N, M)
end

"""
    Circuit

Quantum circuit consisting of a unitary gate chain and measurement operators.
"""
struct Circuit{N}
    moments::Vector{Moment}
    meas::Vector{MeasurementOperator}

    @doc """
        Circuit{N}(mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing)

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `CircuitGate` objects.
    """
    function Circuit{N}(mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing) where {N}
        if isnothing(mops)
            return new{N}(Moment[])
        end
        meas_N = maximum(size.(mops))
        meas_N <= N || error("Measurement operators affecting $meas_N wires provided for Circuit of size $N")
        check_commute(mops) || error("Measurement operators do not commute") 
        return new{N}(Moment[], mops)
    end

    @doc """
        Circuit{N}(gates::Vector{<:CircuitGate}, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing)

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `CircuitGate` objects.
    """
    function Circuit{N}(gates::Vector{<:CircuitGate}, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing) where {N}
        j = 1
        iwires = falses(N)
        moments = Moment[]
        for (i, cg) in enumerate(gates)
            cgwires = collect(cg.iwire)
            all(cgwires .<= N) || error("Unable to add CircuitGate with `iwire`: $(cg.iwire), with maximum circuit size `N`: $N")
            if any(iwires[cgwires])
                push!(moments, Moment(gates[j:i-1]))
                j = i
                iwires[:] .= false
            end
            iwires[cgwires] .= true
        end
        push!(moments, Moment(gates[j:end]))

        if isnothing(mops)
            return new{N}(moments)
        end
        meas_N = maximum(size.(mops))
        meas_N <= N || error("Measurement operators affecting $meas_N wires provided for Circuit of size $N")
        check_commute(mops) || error("Measurement operators do not commute") 
        return new{N}(moments, mops)
    end

    @doc """
        Circuit{N}(gate::AbstractCircuitGate, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing)

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `AbstractCircuitGate` objects.
    """
    function Circuit{N}(gate::CircuitGate, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing) where {N}
        if isnothing(mops)
            return new{N}(Moment(gate))
        end
        new{N}(Moment(gate), mops)
    end

    @doc """
    Circuit{N}(moments::Vector{Moment}, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing) where {N}

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `AbstractCircuitGate` objects.
    """
    function Circuit{N}(moments::Vector{Moment}, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing) where {N}
        if isnothing(mops)
            return new{N}(moments)
        end
        new{N}(moments, mops)
    end
end

matrix(c::Circuit{N}) where {N} = matrix(c.moments, N)

function add_measurement!(c::Circuit{N}, mops::Vector{<:MeasurementOperator}) where {N}
    meas_N = maximum(size.(mops))
    meas_N <= N || error("Measurement operators affecting $meas_N wires provided for Circuit of size $N")
    check_commute(mops) || error("Measurement Operators do not commute") 
    c.meas = mops
end

# make Circuit iterable and indexable
function Base.getindex(c::Circuit, i::Integer) 
    1 <= i <= length(c.moments) || throw(BoundsError(c, i))
    return c.moments[i]
end

function Base.iterate(c::Circuit, state=1) 
    return state > length(c.moments) ? nothing : (c.moments[state], state + 1)
end

# implement methods required for iteration
function Base.firstindex(::Circuit) 
    return 1
end

function Base.lastindex(c::Circuit)
    return length(c.moments)
end

function distribution(c::Circuit, ψ::AbstractVector) 
    return apply(c.moments, ψ)
end

function Base.append!(c::Circuit{N}, gate::CircuitGate) where {N}
    cg_N = maximum(gate.iwire)
    cg_N <= N || error("CircuitGate `cg` has maximum iwire `N`, $cg_N. However, Circuit object `c` has size of $N")
    push!(c.moments, Moment(gate))
end

function Base.append!(c::Circuit{N}, gates::Vector{<:CircuitGate}) where {N}
    j = 1
    iwires = falses(N)
    buffer = CircuitGate[]
    for (i, cg) in enumerate(gates)
        cgwires = collect(cg.iwire)
        all(cgwires .<= N) || error("Unable to add gate with `iwire`: $(cg.iwire), with maximum circuit size `N`: $(N)")
        if any(iwires[cgwires])
            push!(c.moments, Moment(buffer))
            j = i
            iwires[:] .= false
            iwires[cgwires] .= true
            buffer = CircuitGate[cg]
        else
            iwires[cgwires] .= true
            push!(buffer, cg)
        end
    end
    push!(c.moments, Moment(buffer))
end

Base.size(::Circuit{N}) where {N} = N

Base.length(c::Circuit) = length(c.moments)