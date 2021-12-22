"""
    MeasurementOperator{M,G}

General measurement operator. `M` is the number of wires affected by the measurement operator, `G` is the actual operator type (gate or sparse matrix).
```
julia> m = MeasurementOperator(X, (1,))
MeasurementOperator{1,XGate}(XGate(), (1,))

julia> m = MeasurementOperator([0 1; 1 0], (1,))
MeasurementOperator{1,SparseArrays.SparseMatrixCSC{Complex{FloatQ},Int64}}(
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
        new{M,SparseMatrixCSC{Complex{FloatQ},Int}}(sparse(m), iwire)
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

data(m::MeasurementOperator) = data(m.operator)

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


"""
    req_wires(m::MeasurementOperator{M,G})

Minimum number of required qubit wires in a circuit to host the measurement operator `m`.
```jldoctest
julia> meas = mop(X, (4,));
julia> req_wires(meas)
4
```
"""
req_wires(m::MeasurementOperator{M,G}) where {M,G} = maximum(m.iwire)