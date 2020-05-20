"""
    comm(A, B)

Matrix commutator [A, B]. Pre-defined commute methods for specific circuit gate combinations
"""

# Can make more efficient to detect cases when gate wires overlap!
function comm(A::CircuitGate{L, M, G}, B::CircuitGate{N, M, H}) where {L, N, M, G, H}

    cntrlA, gateA = get_controls(A)
    cntrlB, gateB = get_controls(B)

    if length(intersect(gateA, gateB)) == 0
        return true
    end

    if cntrlA == cntrlB && gateA == gateB
        if G == H
            return true
        end
        return comm(A.gate, B.gate)
    end

    mA = Qaintessent.matrix(A)
    mB = Qaintessent.matrix(B)

    if norm(mA*mB - mB*mA) < 1e-13
        return true
    end
    return false

end

function comm(A::AbstractGate{N}, B::AbstractGate{N}) where {N}
    mA = Qaintessent.matrix(A)
    mB = Qaintessent.matrix(B)

    if norm(mA*mB - mB*mA) < 1e-13
        return true
    end
    return false
end

comm(A::AbstractGate{2}, B::AbstractGate{1}) = comm(B::AbstractGate{1}, A::AbstractGate{2})

comm(A::XGate, B::YGate) = false
comm(A::YGate, B::XGate) = false
comm(A::XGate, B::ZGate) = false
comm(A::ZGate, B::XGate) = false
comm(A::YGate, B::ZGate) = false
comm(A::ZGate, B::YGate) = false

comm(A::XGate, B::RxGate) = true
comm(A::RxGate, B::XGate) = true

comm(A::RyGate, B::YGate) = true
comm(A::YGate, B::RyGate) = true

comm(A::RzGate, B::ZGate) = true
comm(A::ZGate, B::RzGate) = true

comm(A::T, B::T) where {T <:AbstractGate} = true
