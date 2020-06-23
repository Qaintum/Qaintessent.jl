"""
    iscommuting(A,B)

Matrix commutator [A,B]. Pre-defined commute methods for specific circuit gate combinations
"""

# Can make more efficient to detect cases when gate wires overlap!
function iscommuting(A::CircuitGate{L,M,G},B::CircuitGate{N,M,H}) where {L,N,M,G,H}

    cntrlA,gateA = get_controls(A)
    cntrlB,gateB = get_controls(B)

    if length(intersect([cntrlA..., gateA...], [cntrlB..., gateB...])) == 0
        return true
    end

    if cntrlA == cntrlB && gateA == gateB
        if G == H
            return true
        end
        return iscommuting(A.gate, B.gate)
    end

    mA = Qaintessent.matrix(A)
    mB = Qaintessent.matrix(B)

    if norm(mA*mB - mB*mA) < 1e-13
        return true
    end
    return false

end

function iscommuting(A::AbstractGate{N},B::AbstractGate{N}) where {N}
    mA = Qaintessent.matrix(A)
    mB = Qaintessent.matrix(B)

    if norm(mA*mB - mB*mA) < 1e-13
        return true
    end
    return false
end

iscommuting(A::AbstractGate{2},B::AbstractGate{1}) = iscommuting(B::AbstractGate{1},A::AbstractGate{2})

iscommuting(A::XGate,B::YGate) = false
iscommuting(A::YGate,B::XGate) = false
iscommuting(A::XGate,B::ZGate) = false
iscommuting(A::ZGate,B::XGate) = false
iscommuting(A::YGate,B::ZGate) = false
iscommuting(A::ZGate,B::YGate) = false

iscommuting(A::XGate,B::RxGate) = true
iscommuting(A::RxGate,B::XGate) = true

iscommuting(A::RyGate,B::YGate) = true
iscommuting(A::YGate,B::RyGate) = true

iscommuting(A::RzGate,B::ZGate) = true
iscommuting(A::ZGate,B::RzGate) = true

iscommuting(A::T,B::T) where {T <:AbstractGate} = true
