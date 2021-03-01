using LinearAlgebra

"""
    iscommuting(A::AbstractMatrix, B::AbstractMatrix)

returns true if [`AbstractMatrix`](@ref) `A` and `B` commute.
"""
function iscommuting(A::AbstractMatrix, B::AbstractMatrix)
    # require compatible dimensions
    if size(A, 2) == size(B, 1) && size(A, 1) == size(B, 2)
        return (A * B ≈ B * A)
    else
        return false
    end
end


# gates of same type do not necessarily commute, see, e.g., MatrixGate

iscommuting(::XGate, ::XGate) = true
iscommuting(::YGate, ::YGate) = true
iscommuting(::ZGate, ::ZGate) = true

iscommuting(::XGate, ::YGate) = false
iscommuting(::YGate, ::XGate) = false
iscommuting(::XGate, ::ZGate) = false
iscommuting(::ZGate, ::XGate) = false
iscommuting(::YGate, ::ZGate) = false
iscommuting(::ZGate, ::YGate) = false

iscommuting(::RxGate, ::RxGate) = true
iscommuting(::RyGate, ::RyGate) = true
iscommuting(::RzGate, ::RzGate) = true

iscommuting(::XGate,  ::RxGate) = true
iscommuting(::RxGate, ::XGate)  = true

iscommuting(::RyGate, ::YGate) = true
iscommuting(::YGate, ::RyGate) = true

iscommuting(::RzGate, ::ZGate) = true
iscommuting(::ZGate, ::RzGate) = true

# entanglement gates commute with each other
iscommuting(::EntanglementXXGate, ::EntanglementXXGate) = true
iscommuting(::EntanglementXXGate, ::EntanglementYYGate) = true
iscommuting(::EntanglementXXGate, ::EntanglementZZGate) = true
iscommuting(::EntanglementYYGate, ::EntanglementXXGate) = true
iscommuting(::EntanglementYYGate, ::EntanglementYYGate) = true
iscommuting(::EntanglementYYGate, ::EntanglementZZGate) = true
iscommuting(::EntanglementZZGate, ::EntanglementXXGate) = true
iscommuting(::EntanglementZZGate, ::EntanglementYYGate) = true
iscommuting(::EntanglementZZGate, ::EntanglementZZGate) = true

# cover MatrixGate, too
iscommuting(A::AbstractGate, B::AbstractGate) = iscommuting(Qaintessent.sparse_matrix(A), Qaintessent.sparse_matrix(B))

"""
    iscommuting(A::CircuitGate{L,ControlledGate{G}}, B::CircuitGate{M,ControlledGate{G}}) where {L,M,G<:AbstractGate}

returns true if controlled [`CircuitGate`](@ref) `A` and `B` commute.
"""
function iscommuting(A::CircuitGate{L,ControlledGate{G}}, B::CircuitGate{M,ControlledGate{G}}) where {L,M,G<:AbstractGate}
    N = maximum((maximum(B.iwire), maximum(A.iwire)))
    # check whether only control wires overlap, if at all
    if (length(intersect(A.iwire[end - L + 1:end], B.iwire)) == 0 &&
        length(intersect(B.iwire[end - M + 1:end], A.iwire)) == 0)
        return true
    end

    if A.iwire == B.iwire
        return iscommuting(A.gate, B.gate)
    end

    # TODO: more fine-grained case switches for partially overlapping 'iwire'
    iscommuting(Qaintessent.sparse_matrix(A, N), Qaintessent.sparse_matrix(B, N))
end


"""
    iscommuting(A::CircuitGate{L,G}, B::CircuitGate{M,H}) where {L,M,G,H}

returns true if general [`CircuitGate`](@ref) `A` and `B` commute.
"""
function iscommuting(A::CircuitGate{L,G}, B::CircuitGate{M,H}) where {L,M,G,H}
    N = maximum((maximum(B.iwire), maximum(A.iwire)))
    if length(intersect(A.iwire, B.iwire)) == 0
        return true
    end

    if A.iwire == B.iwire
        return iscommuting(A.gate, B.gate)
    end

    # TODO: more fine-grained case switches for partially overlapping 'iwire'
    iscommuting(Qaintessent.sparse_matrix(A,N), Qaintessent.sparse_matrix(B,N))
end

"""
    iscommuting(A::MeasurementOperator{L,G}, B::MeasurementOperator{M,H}) where {L,M,G,H}

returns true if [`MeasurementOperator`](@ref) `A` and `B` commute.
"""
function iscommuting(A::MeasurementOperator{L,G}, B::MeasurementOperator{M,H}) where {L,M,G,H}
    N = maximum((maximum(A.iwire), maximum(B.iwire)))
    if length(intersect(A.iwire, B.iwire)) == 0
        return true
    end

    if A.iwire == B.iwire
        return iscommuting(A.operator, B.operator)
    end
    
    iscommuting(sparse_matrix(A, N), sparse_matrix(B, N))
end


function iscommuting(A::MeasurementOperator{M,G}, B::MeasurementOperator{L,H}) where {L,M,G<:AbstractGate,H<:AbstractMatrix}
    if length(intersect(A.iwire, B.iwire)) == 0
        return true
    end

    N = maximum((maximum(A.iwire), maximum(B.iwire)))
    iscommuting(sparse_matrix(A, N), sparse_matrix(B, N))
end

iscommuting(A::MeasurementOperator{M,G}, B::MeasurementOperator{L,H}) where {L,M,G<:AbstractMatrix,H<:AbstractGate} = iscommuting(B, A)