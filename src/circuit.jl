"""
    MeasurementOperator{M,G}

General measurement operator. `M` is the number of wires affected by the measurement operator, `G` is the actual operator type (gate or sparse matrix).
```
julia> m = MeasurementOperator(X, (1,))
MeasurementOperator{1,XGate}(XGate(), (1,))

julia> m = MeasurementOperator([0 1; 1 0], (1,))
MeasurementOperator{1,SparseArrays.SparseMatrixCSC{Complex{Float64},Int64}}(
  [2, 1]  =  1.0+0.0im
  [1, 2]  =  1.0+0.0im, (1,))
```
"""
struct MeasurementOperator{M,G}
    operator::G
    iwire::NTuple{M,Int}
    function MeasurementOperator(g::G, iwire::NTuple{M,Integer}) where {M,G <: AbstractGate}
        M == num_wires(g) || error("CircuitGate affecting $(num_wires(g)) wires given, `iwire` of length $(length(iwire)) provided")
        g == adjoint(g) || error("Measurement operator must be Hermitian.")
        new{M,G}(g, iwire)
    end

    function MeasurementOperator(m::G, iwire::NTuple{M,Integer}) where {M,G <: AbstractMatrix}
        d = 2
        size(m) == (d^M, d^M) || error("Measurement operator must be a $(d^M) × $(d^M) matrix.")
        m ≈ adjoint(m) || error("Measurement operator must be Hermitian.")
        new{M,SparseMatrixCSC{Complex{Float64},Int}}(sparse(m), iwire)
    end

    function MeasurementOperator(m, iwire::Integer...)
        MeasurementOperator(m, iwire)
    end
end

mop(operator, iwire) = MeasurementOperator(operator, iwire)

@memoize function num_wires(m::MeasurementOperator)
    length(m.iwire)
end
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
    sparse_matrix(m::MeasurementOperator{M,G}, N::Integer=0) where {M,G}

returns matrix representation of a `MeasurementOperator{M,G}` object that can applied to a state vector of `N` qubits.
"""
function sparse_matrix(m::MeasurementOperator{M,G}, N::Integer=0) where {M,G}
    # convert to array
    iwire = collect(m.iwire)

    if N == 0
        N = maximum([iwire' M])
    else
        N >= maximum(iwire) || error("CircuitGate applied to iwires, $iwire. Input circuit size `N` is $N")
    end

    gmat = sparse_matrix(m.operator)

    distribute_to_wires(gmat, iwire, N, M)
end


sparse_matrix(g::SparseMatrixCSC{Complex{Float64},Int}) = g


"""
    Circuit{N}

Quantum circuit consisting of a unitary gate chain and measurement operators with `N` qubits.
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
        meas_N = maximum(req_wires.(mops))
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

sparse_matrix(c::Circuit{N}) where {N} = sparse_matrix(c.moments, N)

function add_measurement!(c::Circuit{N}, mops::Vector{<:MeasurementOperator}) where {N}
    meas_N = maximum(req_wires.(mops))
    meas_N <= N || error("Measurement operators affecting $meas_N wires provided for circuit of size $N")
    check_commute(mops) || error("Measurement operators do not commute") 
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

num_wires(::Circuit{N}) where {N} = N

Base.length(c::Circuit) = length(c.moments)

function Base.reverse(c::Circuit)
    c_new = deepcopy(c)
    reverse!(c_new.moments)
    reverse!.(c_new.moments)
    return c_new
end

function Base.reverse!(c::Circuit)
    reverse!(c.moments)
    reverse!.(c.moments)
    return c
end