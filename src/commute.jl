
"""
    iscommuting(A, B)

Test whether two matrices commute.
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

iscommuting(A::XGate, B::XGate) = true
iscommuting(A::YGate, B::YGate) = true
iscommuting(A::ZGate, B::ZGate) = true

iscommuting(A::XGate, B::YGate) = false
iscommuting(A::YGate, B::XGate) = false
iscommuting(A::XGate, B::ZGate) = false
iscommuting(A::ZGate, B::XGate) = false
iscommuting(A::YGate, B::ZGate) = false
iscommuting(A::ZGate, B::YGate) = false

iscommuting(A::RxGate, B::RxGate) = true
iscommuting(A::RyGate, B::RyGate) = true
iscommuting(A::RzGate, B::RzGate) = true

iscommuting(A::XGate,  B::RxGate) = true
iscommuting(A::RxGate, B::XGate)  = true

iscommuting(A::RyGate, B::YGate) = true
iscommuting(A::YGate, B::RyGate) = true

iscommuting(A::RzGate, B::ZGate) = true
iscommuting(A::ZGate, B::RzGate) = true

# entanglement gates commute with each other
iscommuting(A::EntanglementXXGate, B::EntanglementXXGate) = true
iscommuting(A::EntanglementXXGate, B::EntanglementYYGate) = true
iscommuting(A::EntanglementXXGate, B::EntanglementZZGate) = true
iscommuting(A::EntanglementYYGate, B::EntanglementXXGate) = true
iscommuting(A::EntanglementYYGate, B::EntanglementYYGate) = true
iscommuting(A::EntanglementYYGate, B::EntanglementZZGate) = true
iscommuting(A::EntanglementZZGate, B::EntanglementXXGate) = true
iscommuting(A::EntanglementZZGate, B::EntanglementYYGate) = true
iscommuting(A::EntanglementZZGate, B::EntanglementZZGate) = true

# same number of control and target wires
iscommuting(A::ControlledGate{M,N}, B::ControlledGate{M,N}) where {M,N} = iscommuting(A.U, B.U)

# cover MatrixGate, too
iscommuting(A::AbstractGate{N}, B::AbstractGate{N}) where {N} = iscommuting(Qaintessent.matrix(A), Qaintessent.matrix(B))

# catch case of incompatible dimensions
iscommuting(A::AbstractGate{M}, B::AbstractGate{N}) where {M,N} = false


"""
    iscommuting(A, B)

Test whether two controlled circuit gates commute.
"""
function iscommuting(A::CircuitGate{L,N,ControlledGate{S,L}}, B::CircuitGate{M,N,ControlledGate{T,M}}) where {L,M,N,S,T}

    # check whether only control wires overlap, if at all
    if (length(intersect(A.iwire[end - S + 1:end], B.iwire)) == 0 &&
        length(intersect(B.iwire[end - T + 1:end], A.iwire)) == 0)
        return true
    end

    if A.iwire == B.iwire
        return iscommuting(A.gate, B.gate)
    end

    # TODO: more fine-grained case switches for partially overlapping 'iwire'

    iscommuting(Qaintessent.matrix(A), Qaintessent.matrix(B))
end


"""
    iscommuting(A, B)

Test whether two general circuit gates commute.
"""
function iscommuting(A::CircuitGate{L,N,G}, B::CircuitGate{M,N,H}) where {L,M,N,G,H}

    if length(intersect(A.iwire, B.iwire)) == 0
        return true
    end

    if A.iwire == B.iwire
        return iscommuting(A.gate, B.gate)
    end

    # TODO: more fine-grained case switches for partially overlapping 'iwire'

    iscommuting(Qaintessent.matrix(A), Qaintessent.matrix(B))
end


# catch case of different number of wires
iscommuting(A::CircuitGate{K,L,G}, B::CircuitGate{M,N,H}) where {K,L,M,N,G,H} = false
