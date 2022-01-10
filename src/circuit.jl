sparse_matrix(g::SparseMatrixCSC{Complex{FloatQ},Int}) = g


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
    function Circuit{N}(mops::Vector{<:MeasurementOperator}=MeasurementOperator[]) where {N}
        if isempty(mops)
            return new{N}(Moment[], mops)
        end
        mops_N = maximum(req_wires.(mops))
        mops_N <= N || error("Measurement operators affecting $mops_N wires provided for Circuit of size $N")
        check_commute(mops) || error("Measurement operators do not commute") 
        return new{N}(Moment[], mops)
    end

    @doc """
        Circuit{N}(gates::Vector{<:CircuitGate}, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing)

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `CircuitGate` objects.
    """
    function Circuit{N}(gates::Vector{<:CircuitGate}, mops::Vector{<:MeasurementOperator}=MeasurementOperator[]) where {N}
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
        
        if isempty(mops)
            return new{N}(moments, mops)
        end
        mops_N = maximum(size.(mops))
        mops_N <= N || error("Measurement operators affecting $mops_N wires provided for Circuit of size $N")
        check_commute(mops) || error("Measurement operators do not commute") 
        return new{N}(moments, mops)
    end

    @doc """
        Circuit{N}(gate::AbstractCircuitGate, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing)

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `AbstractCircuitGate` objects.
    """
    function Circuit{N}(gate::CircuitGate, mops::Vector{<:MeasurementOperator}=MeasurementOperator[]) where {N}
        if isempty(mops)
            return new{N}([Moment(gate)], mops)
        end

        mops_N = maximum(size.(mops))
        mops_N <= N || error("Measurement operators affecting $mops_N wires provided for Circuit of size $N")
        check_commute(mops) || error("Measurement operators do not commute") 
        new{N}([Moment(gate)], mops)
    end

    @doc """
    Circuit{N}(moments::Vector{Moment}, mops::Union{Vector{<:MeasurementOperator},Nothing}=nothing) where {N}

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `AbstractCircuitGate` objects.
    """
    function Circuit{N}(moments::Vector{Moment}, mops::Vector{<:MeasurementOperator}=MeasurementOperator[]) where {N}
        if isempty(mops)
            return new{N}(moments, mops)
        end
        
        mops_N = maximum(size.(mops))
        mops_N <= N || error("Measurement operators affecting $mops_N wires provided for Circuit of size $N")
        check_commute(mops) || error("Measurement operators do not commute") 

        new{N}(moments, mops)
    end
end

sparse_matrix(c::Circuit{N}) where {N} = sparse_matrix(c.moments, N)

function add_measurement!(c::Circuit{N}, mops::Vector{<:MeasurementOperator}) where {N}
    if isempty(mops)
        return c
    end
    meas_N = maximum(req_wires.(mops))
    meas_N <= N || error("Measurement operators affecting $meas_N wires provided for circuit of size $N")
    
    check_commute([c.meas; mops]) || error("Measurement operators do not commute") 
    
    append!(c.meas, mops)
    return c
end

function add_measurement!(c::Circuit{N}, mop::MeasurementOperator) where {N}
    meas_N = maximum(req_wires(mop))
    meas_N <= N || error("Measurement operators affecting $meas_N wires provided for circuit of size $N")
    check_commute([c.meas; mop]) || error("Measurement operators do not commute") 
    
    push!(c.meas, mop)
end

# make Circuit iterable and indexable
function Base.getindex(c::Circuit, i::Integer) 
    1 <= i <= length(c.moments) || throw(BoundsError(c, i))
    return c.moments[i]
end

function Base.getindex(c::Circuit, ::Colon) 
    return c.moments[:]
end

function Base.getindex(c::Circuit, ur::UnitRange{Int64}) 
    return c.moments[ur]
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

function distribution(ψ::AbstractVector, c::Circuit) 
    return apply(ψ, c.moments)
end

function Base.append!(c::Circuit{N}, gate::CircuitGate) where {N}
    cg_N = maximum(gate.iwire)
    cg_N <= N || error("CircuitGate `cg` has maximum iwire `N`, $cg_N. However, Circuit object `c` has size of $N")
    push!(c.moments, Moment(gate))
end

function Base.append!(c::Circuit{N}, gates::Vector{<:CircuitGate}) where {N}
    j = 1
    iwires = falses(N)
    if !isempty(c)
        buffer = pop!(c).gates

        for cg in buffer
            cgwires = collect(cg.iwire)
            iwires[cgwires] .= true
        end
    else
        buffer = CircuitGate[]
    end

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

function Base.pop!(c::Circuit)
    pop!(c.moments)
end

function Base.isempty(c::Circuit)
    isempty(c.moments)
end