
"""Tailored conjugation of density matrix by XGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,XGate})
    # qubit index the gate acts on
    j = 2*cg.iwire[1]-1
    v = deepcopy(ρ.v)
    for i in 0:4^ρ.N-1
        if (i >> j) & 1 == 1
            @inbounds v[i+1] = -ρ.v[i+1]
        end
    end
    
    return DensityMatrix(v, ρ.N)
end

"""Tailored conjugation of density matrix by XGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,XGate})
    # qubit index the gate acts on
    j = 2*cg.iwire[1]-1
    for i in 0:4^ρ.N-1
        if (i >> j) & 1 == 1
            @inbounds ρ.v[i+1] = -ρ.v[i+1]
        end
    end
    
    return ρ
end

"""Tailored implementation of 1/2 (X ρ + ρ X)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,XGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρ.scratch .= ρ.v
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    vs = reshape(ρ.scratch, 4^(j-1), 4, 4^(ρ.N-j))

    vs[:, 1, :] .= ρv[:, 2, :]      # 1/2 (X X + X X) = I
    vs[:, 2, :] .= ρv[:, 1, :]      # 1/2 (X I + I X) = X
    vs[:, 3, :] .= 0                # 1/2 (X Y + Y X) = 0
    vs[:, 4, :] .= 0                # 1/2 (X Z + Z X) = 0

    return DensityMatrix(reshape(vs,:), ρ.N)
end

"""Tailored implementation of 1/2 (X ρ + ρ X)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,XGate}, factor::FloatQ=convert(FloatQ, 0))
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    @inbounds ρ.scratch .= factor .* ρ.v
    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 0
            @inbounds ρ.scratch[i+1] += ρ.v[i+1+shift]
        elseif (i >> l) & 3 == 1
            @inbounds ρ.scratch[i+1] += ρ.v[i+1-shift]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (X ρ - ρ X)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,XGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  0               # i/2 (X I - I X) =  0
    vs[:, 2, :] .=  0               # i/2 (X X - X X) =  0
    vs[:, 3, :] .=  ρv[:, 4, :]     # i/2 (X Z - Z X) =  Y
    vs[:, 4, :] .= -ρv[:, 3, :]     # i/2 (X Y - Y X) = -Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (X ρ - ρ X)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,XGate})
    # qubit index the gate acts on
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    shift = 2^l

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            @inbounds ρ.scratch[i+1] = 0
        elseif (i >> l) & 1 == 0
            @inbounds ρ.scratch[i+1] = ρ.v[i+shift+1]
        else
            @inbounds ρ.scratch[i+1] = -ρ.v[i-shift+1]
        end
    end
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by YGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,YGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # Y I Y =  I
    vs[:, 2, :] .= -ρv[:, 2, :]     # Y X Y = -X
    vs[:, 3, :] .=  ρv[:, 3, :]     # Y Y Y =  Y
    vs[:, 4, :] .= -ρv[:, 4, :]     # Y Z Y = -Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by YGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,YGate})
    # qubit index the gate acts on
    l = 2*cg.iwire[1]-2
    shift = 2^l
    for i in 0:4^ρ.N-1
        if (i >> l) & 1 == 1
            @inbounds ρ.v[i+1] = -ρ.v[i+1]
        end
    end    

    return ρ
end

"""Tailored implementation of 1/2 (Y ρ + ρ Y)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,YGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .= ρv[:, 3, :]      # 1/2 (Y Y + Y Y) = I
    vs[:, 2, :] .= 0                # 1/2 (Y X + X Y) = 0
    vs[:, 3, :] .= ρv[:, 1, :]      # 1/2 (Y I + I Y) = Y
    vs[:, 4, :] .= 0                # 1/2 (Y Z + Z Y) = 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Y ρ + ρ Y)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,YGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^k
    @inbounds ρ.scratch .= factor * ρ.v
    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 0
            @inbounds ρ.scratch[i+1] += ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 2
            @inbounds ρ.scratch[i+1] += ρ.v[i-shift+1]
        end
    end    
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (Y ρ - ρ Y)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,YGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  0               # i/2 (Y I - I Y) =  0
    vs[:, 2, :] .= -ρv[:, 4, :]     # i/2 (Y Z - Z Y) = -X
    vs[:, 3, :] .=  0               # i/2 (Y Y - Y Y) =  0
    vs[:, 4, :] .=  ρv[:, 2, :]     # i/2 (Y X - X Y) =  Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Y ρ - ρ Y)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,YGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^k
    for i in 0:4^ρ.N-1
        if (i >> l) & 1 == 0
            @inbounds ρ.scratch[i+1] = 0
        elseif (i >> k) & 1 == 1
            @inbounds ρ.scratch[i+1] = ρ.v[i-shift+1]
        else
            @inbounds ρ.scratch[i+1] = -ρ.v[i+shift+1]
        end
    end    
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by ZGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,ZGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # Z I Z =  I
    vs[:, 2, :] .= -ρv[:, 2, :]     # Z X Z = -X
    vs[:, 3, :] .= -ρv[:, 3, :]     # Z Y Z = -Y
    vs[:, 4, :] .=  ρv[:, 4, :]     # Z Z Z =  Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by ZGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,ZGate})
    # qubit index the gate acts on
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 2 || (i >> l) & 3 == 1
            @inbounds ρ.v[i+1] = -ρ.v[i+1]
        end
    end    

    return ρ
end

"""Tailored implementation of 1/2 (Z ρ + ρ Z)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,ZGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .= ρv[:, 4, :]      # 1/2 (Z Z + Z Z) = I
    vs[:, 2, :] .= 0                # 1/2 (Z X + X Z) = 0
    vs[:, 3, :] .= 0                # 1/2 (Z Y + Y Z) = 0
    vs[:, 4, :] .= ρv[:, 1, :]      # 1/2 (Z I + I Z) = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Z ρ + ρ Z)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,ZGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^k + 2^l

    @inbounds ρ.scratch .= factor * ρ.v

    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 0
            @inbounds ρ.scratch[i+1] += ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 3
            @inbounds ρ.scratch[i+1] += ρ.v[i-shift+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Z ρ - ρ Z)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,ZGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  0               # i/2 (Z I - I Z) =  0
    vs[:, 2, :] .=  ρv[:, 3, :]     # i/2 (Z Y - Y Z) =  X
    vs[:, 3, :] .= -ρv[:, 2, :]     # i/2 (Z X - X Z) = -Y
    vs[:, 4, :] .=  0               # i/2 (Z Z - Z Z) =  0

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Z ρ - ρ Z)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,ZGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    for i in 0:4^ρ.N-1
        if (i >> l) & 1 == (i >> k) & 1
            @inbounds ρ.scratch[i+1] = 0
        elseif (i >> l) & 1 == 1
            @inbounds ρ.scratch[i+1] = ρ.v[i+shift+1]
        else
            @inbounds ρ.scratch[i+1] = -ρ.v[i-shift+1]
        end
    end    

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by the Hadamard gate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,HadamardGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # H I H =  I
    vs[:, 2, :] .=  ρv[:, 4, :]     # H Z H =  X
    vs[:, 3, :] .= -ρv[:, 3, :]     # H Y H = -Y
    vs[:, 4, :] .=  ρv[:, 2, :]     # H X H =  Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by the Hadamard gate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,HadamardGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^k

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = ρ.v[i+1]
            else
                @inbounds ρ.scratch[i+1] = ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = -ρ.v[i+1]
            else
                @inbounds ρ.scratch[i+1] = ρ.v[i-shift+1]
            end
        end
    end
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of 1/2 (H ρ + ρ H)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,HadamardGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    # 1/2 (H X + X H) = I/√2
    # 1/2 (H Z + Z H) = I/√2
    # 1/2 (H I + I H) = H = (X + Z)/√2
    vs = similar(ρv)
    vs[:, 1, :] .= (ρv[:, 2, :] .+ ρv[:, 4, :]) ./ sqrt(2)
    vs[:, 2, :] .= ρv[:, 1, :] ./ sqrt(2)
    vs[:, 3, :] .= 0
    vs[:, 4, :] .= ρv[:, 1, :] ./ sqrt(2)

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (H ρ + ρ H)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,HadamardGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    @inbounds ρ.scratch .= factor .* ρ.v
    for i in 0:4^ρ.N-1
        if (i >> l) & 1 == 1
            bit = (i >> k) & 1
            @inbounds ρ.scratch[i+1] += ρ.v[i-(1+2bit)*shift+1] / sqrt(2)
        elseif (i >> k) & 1 == 0
            @inbounds ρ.scratch[i+1] += (ρ.v[i+3shift+1] + ρ.v[i+shift+1]) / sqrt(2)
        end
    end
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (H ρ - ρ H)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,HadamardGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    # i/2 (H X - X H) = -Y/√2
    # i/2 (H Y - Y H) = (X - Z)/√2
    # i/2 (H Z - Z H) =  Y/√2
    vs = similar(ρv)
    vs[:, 1, :] .=  0
    vs[:, 2, :] .=  ρv[:, 3, :] ./ sqrt(2)
    vs[:, 3, :] .= (ρv[:, 4, :] .- ρv[:, 2, :]) ./ sqrt(2)
    vs[:, 4, :] .= -ρv[:, 3, :] ./ sqrt(2)

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (H ρ - ρ H)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,HadamardGate})
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = 0
            else
                @inbounds ρ.scratch[i+1] = ρ.v[i+shift+1] / sqrt(2)
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = (ρ.v[i+shift+1] - ρ.v[i-shift+1]) / sqrt(2)
            else
                @inbounds ρ.scratch[i+1] = -ρ.v[i-shift+1] / sqrt(2)
            end
        end
    end
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by SGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,SGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # S I S† =  I
    vs[:, 2, :] .= -ρv[:, 3, :]     # S Y S† = -X
    vs[:, 3, :] .=  ρv[:, 2, :]     # S X S† =  Y
    vs[:, 4, :] .=  ρv[:, 4, :]     # S Z S† =  Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by SGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,SGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 1
            @inbounds ρ.scratch[i+1] = -ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 2
            @inbounds ρ.scratch[i+1] = ρ.v[i-shift+1]
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end    

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of 1/2 (S ρ + ρ S†)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,SGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = 0.5 .* ρv
    vs[:, 1, :] .+= 0.5 .* ρv[:, 4, :]
    vs[:, 2, :] .-= 0.5 .* ρv[:, 3, :]
    vs[:, 3, :] .+= 0.5 .* ρv[:, 2, :]
    vs[:, 4, :] .+= 0.5 .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (S ρ + ρ S†)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,SGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    factor = factor + 0.5
    @inbounds ρ.scratch .= factor .* ρ.v
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] -= 0.5ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (S ρ - ρ S†)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,SGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = -0.5 .* ρv
    vs[:, 1, :] .+= 0.5 .* ρv[:, 4, :]
    vs[:, 2, :] .+= 0.5 .* ρv[:, 3, :]
    vs[:, 3, :] .-= 0.5 .* ρv[:, 2, :]
    vs[:, 4, :] .+= 0.5 .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (S ρ - ρ S†)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,SGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    @inbounds ρ.scratch .= -0.5 .* ρ.v
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= 0.5ρ.v[i-shift+1]
            else        
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored conjugation of density matrix by SdagGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,SdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # S† I S =  I
    vs[:, 2, :] .=  ρv[:, 3, :]     # S† Y S =  X
    vs[:, 3, :] .= -ρv[:, 2, :]     # S† X S = -Y
    vs[:, 4, :] .=  ρv[:, 4, :]     # S† Z S =  Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by SdagGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,SdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 1
            @inbounds ρ.scratch[i+1] = ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 2
                @inbounds ρ.scratch[i+1] = -ρ.v[i-shift+1]
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (S† ρ + ρ S)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,SdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = 0.5 .* ρv
    vs[:, 1, :] .+= 0.5 .* ρv[:, 4, :]
    vs[:, 2, :] .+= 0.5 .* ρv[:, 3, :]
    vs[:, 3, :] .-= 0.5 .* ρv[:, 2, :]
    vs[:, 4, :] .+= 0.5 .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (S† ρ + ρ S)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,SdagGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    factor = factor + 0.5
    @inbounds ρ.scratch .= factor .* ρ.v
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= 0.5ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += 0.5ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (S† ρ - ρ S)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,SdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = 0.5 .* ρv
    vs[:, 1, :] .-= 0.5 .* ρv[:, 4, :]
    vs[:, 2, :] .+= 0.5 .* ρv[:, 3, :]
    vs[:, 3, :] .-= 0.5 .* ρv[:, 2, :]
    vs[:, 4, :] .-= 0.5 .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (S† ρ - ρ S)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,SdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = -0.5ρ.v[i+3shift+1] + 0.5ρ.v[i+1]
            else
                @inbounds ρ.scratch[i+1] = 0.5ρ.v[i+shift+1] + 0.5ρ.v[i+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = -0.5ρ.v[i-shift+1] + 0.5ρ.v[i+1]
            else
                @inbounds ρ.scratch[i+1] = -0.5ρ.v[i-3shift+1] + 0.5ρ.v[i+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by TGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,TGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]                             # T I T† = I
    vs[:, 2, :] .= (ρv[:, 2, :] .- ρv[:, 3, :]) ./ sqrt(2)  # (T X T† - T Y T†)/√2 = X
    vs[:, 3, :] .= (ρv[:, 2, :] .+ ρv[:, 3, :]) ./ sqrt(2)  # (T X T† + T Y T†)/√2 = Y
    vs[:, 4, :] .=  ρv[:, 4, :]                             # T Z T† = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by TGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,TGate})
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 1
            @inbounds ρ.scratch[i+1] = (ρ.v[i+1] - ρ.v[i+shift+1]) / sqrt(2)
        elseif (i >> l) & 3 == 2
            @inbounds ρ.scratch[i+1] = (ρ.v[i+1] + ρ.v[i-shift+1]) / sqrt(2)
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of 1/2 (T ρ + ρ T†)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,TGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = (0.5 * (1 + 1/sqrt(2))) .* ρv
    vs[:, 1, :] .+= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 4, :]
    vs[:, 2, :] .-= (0.5/sqrt(2))           .* ρv[:, 3, :]
    vs[:, 3, :] .+= (0.5/sqrt(2))           .* ρv[:, 2, :]
    vs[:, 4, :] .+= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (T ρ + ρ T†)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,TGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    n1 = (0.5 * (1 + 1/sqrt(2))) + factor
    n2 = (0.5 * (1 - 1/sqrt(2)))
    n3 = 0.5/sqrt(2)

    @inbounds ρ.scratch .= n1 .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += n2 * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] -= n3 * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += n3 * ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += n2 * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (T ρ - ρ T†)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,TGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = (-0.5/sqrt(2)) .* ρv
    vs[:, 1, :] .+= (0.5/sqrt(2))           .* ρv[:, 4, :]
    vs[:, 2, :] .+= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 3, :]
    vs[:, 3, :] .-= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 2, :]
    vs[:, 4, :] .+= (0.5/sqrt(2))           .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (T ρ - ρ T†)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,TGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    n1 = -0.5/sqrt(2)
    n2 = 0.5/sqrt(2)
    n3 = 0.5 * (1 - 1/sqrt(2))

    @inbounds ρ.scratch .= n1 .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += n2 * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] += n3 * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= n3 * ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += n2 * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Tailored conjugation of density matrix by TdagGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,TdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=   ρv[:, 1, :]                            # T† I T = I
    vs[:, 2, :] .= ( ρv[:, 2, :] .+ ρv[:, 3, :]) ./ sqrt(2) # ( T† X T + T† Y T)/√2 = X
    vs[:, 3, :] .= (-ρv[:, 2, :] .+ ρv[:, 3, :]) ./ sqrt(2) # (-T† X T + T† Y T)/√2 = Y
    vs[:, 4, :] .=   ρv[:, 4, :]                            # T† Z T = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by TdagGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,TdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l
    
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 != (i >> l) & 1
            if (i >> l) & 1 == 1
                @inbounds ρ.scratch[i+1] = (ρ.v[i+1] + ρ.v[i+shift+1]) / sqrt(2)
            else
                @inbounds ρ.scratch[i+1] = (ρ.v[i+1] - ρ.v[i-shift+1]) / sqrt(2)
            end
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (T† ρ + ρ T)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,TdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = (0.5 * (1 + 1/sqrt(2))) .* ρv
    vs[:, 1, :] .+= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 4, :]
    vs[:, 2, :] .+= (0.5/sqrt(2))           .* ρv[:, 3, :]
    vs[:, 3, :] .-= (0.5/sqrt(2))           .* ρv[:, 2, :]
    vs[:, 4, :] .+= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (T† ρ + ρ T)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,TdagGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    n1 = 0.5 * (1 + 1/sqrt(2)) + factor
    n2 = 0.5 * (1 - 1/sqrt(2))
    n3 = 0.5/sqrt(2)

    @inbounds ρ.scratch .= n1 .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += n2 * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] += n3 * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= n3 * ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += n2 * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (T† ρ - ρ T)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,TdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = (0.5/sqrt(2)) .* ρv
    vs[:, 1, :] .-= (0.5/sqrt(2))           .* ρv[:, 4, :]
    vs[:, 2, :] .+= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 3, :]
    vs[:, 3, :] .-= (0.5 * (1 - 1/sqrt(2))) .* ρv[:, 2, :]
    vs[:, 4, :] .-= (0.5/sqrt(2))           .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (T† ρ - ρ T)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,TdagGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    n1 = 0.5/sqrt(2)
    n2 = 0.5/sqrt(2)
    n3 = 0.5 * (1 - 1/sqrt(2))

    @inbounds ρ.scratch .= n1 .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= n2 * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] += n3 * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= n3 * ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] -= n2 * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by RxGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,RxGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # Rx(θ) I Rx(-θ) = I
    vs[:, 2, :] .=       ρv[:, 2, :]                        # Rx(θ) X Rx(-θ) = X
    vs[:, 3, :] .= cosθ.*ρv[:, 3, :] .- sinθ.*ρv[:, 4, :]   # cos(θ) Rx(θ) Y Rx(-θ) - sin(θ) Rx(θ) Z Rx(-θ) = Y
    vs[:, 4, :] .= sinθ.*ρv[:, 3, :] .+ cosθ.*ρv[:, 4, :]   # sin(θ) Rx(θ) Y Rx(-θ) + cos(θ) Rx(θ) Z Rx(-θ) = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end


"""Tailored conjugation of density matrix by RxGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,RxGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])
    
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 1
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] - sinθ * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] + sinθ * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (Rx(θ) ρ + ρ Rx(-θ))"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,RxGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = cosθ2 .* ρv
    vs[:, 3, :] .-= sinθ2 .* ρv[:, 4, :]
    vs[:, 4, :] .+= sinθ2 .* ρv[:, 3, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Rx(θ) ρ + ρ Rx(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,RxGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v
    
    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 2
            @inbounds ρ.scratch[i+1] -= sinθ2 * ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 3
            @inbounds ρ.scratch[i+1] += sinθ2 * ρ.v[i-shift+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Rx(θ) ρ - ρ Rx(-θ))"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,RxGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = similar(ρv)
    vs[:, [1, 2], :] .= sinθ2 .* ρv[:, [2, 1], :]
    vs[:, [3, 4], :] .= 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Rx(θ) ρ - ρ Rx(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,RxGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    sinθ2 = sin(0.5*cg.gate.θ[])
    
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = sinθ2 * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = sinθ2 * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored conjugation of density matrix by RyGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,RyGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # Ry(θ) I Ry(-θ) = I
    vs[:, 2, :] .= sinθ.*ρv[:, 4, :] .+ cosθ.*ρv[:, 2, :]   # sin(θ) Ry(θ) Z Ry(-θ) + cos(θ) Ry(θ) X Ry(-θ) = X
    vs[:, 3, :] .=       ρv[:, 3, :]                        # Ry(θ) Y Ry(-θ) = Y
    vs[:, 4, :] .= cosθ.*ρv[:, 4, :] .- sinθ.*ρv[:, 2, :]   # cos(θ) Ry(θ) Z Ry(-θ) - sin(θ) Ry(θ) X Ry(-θ) = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by RyGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,RyGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^k

    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])
    
    for i in 0:4^ρ.N-1
        if (i >> l) & 1 == 1
            if (i >> k) & 1 == 0
                @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] + sinθ * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] - sinθ * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (Ry(θ) ρ + ρ Ry(-θ))"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,RyGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = cosθ2 .* ρv
    vs[:, 2, :] .+= sinθ2 .* ρv[:, 4, :]
    vs[:, 4, :] .-= sinθ2 .* ρv[:, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Ry(θ) ρ + ρ Ry(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,RyGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    shift = 2^k

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v    
    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 1
            @inbounds ρ.scratch[i+1] += sinθ2 * ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 3
            @inbounds ρ.scratch[i+1] -= sinθ2 * ρ.v[i-shift+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Ry(θ) ρ - ρ Ry(-θ))"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,RyGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = similar(ρv)
    vs[:, [1, 3], :] .= sinθ2 .* ρv[:, [3, 1], :]
    vs[:, [2, 4], :] .= 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Ry(θ) ρ - ρ Ry(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,RyGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    sinθ2 = sin(0.5*cg.gate.θ[])

    shift = 2^k

    for i in 0:4^ρ.N-1
        if (i >> l) & 1 == 0
            if (i >> k) & 1 == 0
                @inbounds ρ.scratch[i+1] = sinθ2 * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = sinθ2 * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Tailored conjugation of density matrix by RzGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,RzGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # Rz(θ) I Rz(-θ) = I
    vs[:, 2, :] .= cosθ.*ρv[:, 2, :] .- sinθ.*ρv[:, 3, :]   # cos(θ) Rz(θ) X Rz(-θ) - sin(θ) Rz(θ) Y Rz(-θ) = X
    vs[:, 3, :] .= sinθ.*ρv[:, 2, :] .+ cosθ.*ρv[:, 3, :]   # sin(θ) Rz(θ) X Rz(-θ) + cos(θ) Rz(θ) Y Rz(-θ) = Y
    vs[:, 4, :] .=       ρv[:, 4, :]                        # Rz(θ) Z Rz(-θ) = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by RzGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,RzGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    shift = 2^l    

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 != (i >> l) & 1
            if (i >> l) & 1 == 1
                @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] -  sinθ * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] +  sinθ * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (Rz(θ) ρ + ρ Rz(-θ))"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,RzGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = cosθ2 .* ρv
    vs[:, 2, :] .-= sinθ2 .* ρv[:, 3, :]
    vs[:, 3, :] .+= sinθ2 .* ρv[:, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Rz(θ) ρ + ρ Rz(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,RzGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    shift = 2^l    

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> l) & 3 == 1
            @inbounds ρ.scratch[i+1] -= sinθ2 * ρ.v[i+shift+1]
        elseif (i >> l) & 3 == 2
            @inbounds ρ.scratch[i+1] += sinθ2 * ρ.v[i-shift+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Rz(θ) ρ - ρ Rz(-θ))"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,RzGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = similar(ρv)
    vs[:, [1, 4], :] .= sinθ2 .* ρv[:, [4, 1], :]
    vs[:, [2, 3], :] .= 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Rz(θ) ρ - ρ Rz(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,RzGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    sinθ2 = sin(0.5*cg.gate.θ[])

    shift = 3*2^l    

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == (i >> l) & 1
            if (i >> k) & 1 == 0
                @inbounds ρ.scratch[i+1] = sinθ2 * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = sinθ2 * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored conjugation of density matrix by RotationGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,RotationGate})
    θ = norm(cg.gate.nθ)
    if θ == 0
        # for consistency, we return a copy here
        return DensityMatrix(copy(ρ.v), ρ.N)
    end
    cosθ = cos(θ)
    sinθ = sin(θ)
    n = cg.gate.nθ/θ

    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .= ρv[:, 1, :]      # Rn(θ) I Rn(-θ) = I
    # Rodrigues' rotation formula
    vs[:, 2, :] .= cosθ.*ρv[:, 2, :] .+ sinθ.*(n[2].*ρv[:, 4, :] .- n[3].*ρv[:, 3, :]) .+ ((1 - cosθ)*n[1]).*(n[1].*ρv[:, 2, :] .+ n[2].*ρv[:, 3, :] .+ n[3].*ρv[:, 4, :])
    vs[:, 3, :] .= cosθ.*ρv[:, 3, :] .+ sinθ.*(n[3].*ρv[:, 2, :] .- n[1].*ρv[:, 4, :]) .+ ((1 - cosθ)*n[2]).*(n[1].*ρv[:, 2, :] .+ n[2].*ρv[:, 3, :] .+ n[3].*ρv[:, 4, :])
    vs[:, 4, :] .= cosθ.*ρv[:, 4, :] .+ sinθ.*(n[1].*ρv[:, 3, :] .- n[2].*ρv[:, 2, :]) .+ ((1 - cosθ)*n[3]).*(n[1].*ρv[:, 2, :] .+ n[2].*ρv[:, 3, :] .+ n[3].*ρv[:, 4, :])

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by RotationGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,RotationGate})
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    shift = 2^l

    θ = norm(cg.gate.nθ)
    if θ == 0
        # for consistency, we return a copy here
        return ρ
    end
    cosθ = cos(θ)
    sinθ = sin(θ)
    n = cg.gate.nθ/θ

    n21 = cosθ + (1-cosθ)*n[1]^2
    n22 = -n[3]*sinθ + (1-cosθ)*n[1]*n[2]
    n23 = n[2]*sinθ + (1-cosθ)*n[1]*n[3]

    n31 = n[3]*sinθ + (1-cosθ)*n[1]*n[2]
    n32 = cosθ + (1-cosθ)*n[2]^2
    n33 = -n[1]*sinθ + (1-cosθ)*n[2]*n[3]

    n41 = -n[2]*sinθ + (1-cosθ)*n[1]*n[3]
    n42 = n[1]*sinθ + (1-cosθ)*n[2]*n[3]
    n43 = cosθ + (1-cosθ)*n[3]^2

    # qubit index the gate acts on
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = ρ.v[i+1]
            else
                @inbounds ρ.scratch[i+1] = n21 * ρ.v[i+1] + n22 * ρ.v[i+shift+1] + n23 * ρ.v[i+2shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = n31 * ρ.v[i-shift+1] + n32 * ρ.v[i+1] + n33 * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = n41 * ρ.v[i-2shift+1] + n42 * ρ.v[i-shift+1] + n43 * ρ.v[i+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of 1/2 (Rn(θ) ρ + ρ Rn(-θ))"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,RotationGate})
    θ = norm(cg.gate.nθ)
    if θ == 0
        # for consistency, we return a copy here
        return DensityMatrix(copy(ρ.v), ρ.N)
    end
    n = cg.gate.nθ/θ

    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    
    cosθ2 = cos(0.5*θ)
    sinθ2 = sin(0.5*θ)
    n = cg.gate.nθ/θ
    sn = sinθ2 * n

    vs = cosθ2 .* ρv
    vs[:, 2, :] .+= (sn[2] .* ρv[:, 4, :] .- sn[3] .* ρv[:, 3, :])
    vs[:, 3, :] .+= (sn[3] .* ρv[:, 2, :] .- sn[1] .* ρv[:, 4, :])
    vs[:, 4, :] .+= (sn[1] .* ρv[:, 3, :] .- sn[2] .* ρv[:, 2, :])

    return DensityMatrix(reshape(vs, :), ρ.N)
end


"""Tailored implementation of 1/2 (Rn(θ) ρ + ρ Rn(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,RotationGate}, factor::FloatQ=convert(FloatQ, 0))
    θ = norm(cg.gate.nθ)
    if θ == 0
        # for consistency, we return a copy here
        return ρ
    end
    n = cg.gate.nθ/θ
    cosθ2 = cos(0.5*θ)
    sinθ2 = sin(0.5*θ)
    sn = sinθ2 * n

    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    shift = 2^l

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v
    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 1
                @inbounds ρ.scratch[i+1] += sn[2] * ρ.v[i+2shift+1] - sn[3] * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += sn[3] * ρ.v[i-shift+1] - sn[1] * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] += sn[1] * ρ.v[i-shift+1] - sn[2] * ρ.v[i-2shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (Rn(θ) ρ - ρ Rn(-θ))"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,RotationGate})
    θ = norm(cg.gate.nθ)
    if θ == 0
        return DensityMatrix(zero(ρ.v), ρ.N)
    end
    n = cg.gate.nθ/θ

    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    sn = sin(0.5*θ) * n

    vs = similar(ρv)
    vs[:, 1, :] .= sn[1] .* ρv[:, 2, :] + sn[2] .* ρv[:, 3, :] + sn[3] .* ρv[:, 4, :]
    vs[:, 2, :] .= sn[1] .* ρv[:, 1, :]
    vs[:, 3, :] .= sn[2] .* ρv[:, 1, :]
    vs[:, 4, :] .= sn[3] .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Rn(θ) ρ - ρ Rn(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,RotationGate})
    θ = norm(cg.gate.nθ)
    if θ == 0
        return ρ
    end
    n = cg.gate.nθ/θ
    sn = sin(0.5*θ) * n

    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    shift = 2^l

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = sn[1] * ρ.v[i+shift+1] + sn[2] * ρ.v[i+2shift+1] + sn[3] * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] = sn[1] * ρ.v[i-shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] = sn[2] * ρ.v[i-2shift+1]
            else
                @inbounds ρ.scratch[i+1] = sn[3] * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by PhaseShiftGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,PhaseShiftGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosϕ = cos(cg.gate.ϕ[])
    sinϕ = sin(cg.gate.ϕ[])

    # agrees with Rz(θ) since global prefactor cancels
    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # P(ϕ) I P(-ϕ) = I
    vs[:, 2, :] .= cosϕ.*ρv[:, 2, :] .- sinϕ.*ρv[:, 3, :]   # cos(ϕ) P(ϕ) X P(-ϕ) - sin(ϕ) P(ϕ) Y P(-ϕ) = X
    vs[:, 3, :] .= sinϕ.*ρv[:, 2, :] .+ cosϕ.*ρv[:, 3, :]   # sin(ϕ) P(ϕ) X P(-ϕ) + cos(ϕ) P(ϕ) Y P(-ϕ) = Y
    vs[:, 4, :] .=       ρv[:, 4, :]                        # P(ϕ) Z P(-ϕ) = Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by PhaseShiftGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,PhaseShiftGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    cosϕ = cos(cg.gate.ϕ[])
    sinϕ = sin(cg.gate.ϕ[])

    shift = 2^l

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 != (i >> l) & 1
            if (i >> l) & 1 == 1
                @inbounds ρ.scratch[i+1] = cosϕ * ρ.v[i+1] - sinϕ * ρ.v[i+shift+1]
            else
                @inbounds ρ.scratch[i+1] = cosϕ * ρ.v[i+1] + sinϕ * ρ.v[i-shift+1]
            end
        else
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of 1/2 (P(ϕ) ρ + ρ P(-ϕ))"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,PhaseShiftGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    cosϕ2 = cos(0.5*cg.gate.ϕ[])
    sinϕ2 = sin(0.5*cg.gate.ϕ[])

    vs = (cosϕ2^2) .* ρv
    vs[:, 1, :] .+= (sinϕ2^2)     .* ρv[:, 4, :]
    vs[:, 2, :] .-= (sinϕ2*cosϕ2) .* ρv[:, 3, :]
    vs[:, 3, :] .+= (sinϕ2*cosϕ2) .* ρv[:, 2, :]
    vs[:, 4, :] .+= (sinϕ2^2)     .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (P(ϕ) ρ + ρ P(-ϕ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,PhaseShiftGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    sinϕ2cosϕ2 = cos(0.5*cg.gate.ϕ[]) * sin(0.5*cg.gate.ϕ[])
    sinϕ2sinϕ2 = sin(0.5*cg.gate.ϕ[]) * sin(0.5*cg.gate.ϕ[])
    cosϕ2cosϕ2 = cos(0.5*cg.gate.ϕ[]) * cos(0.5*cg.gate.ϕ[])

    shift = 2^l
    @inbounds ρ.scratch .= (cosϕ2cosϕ2+factor) .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += sinϕ2sinϕ2 * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] -= sinϕ2cosϕ2 * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += sinϕ2cosϕ2 * ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += sinϕ2sinϕ2 * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (P(ϕ) ρ - ρ P(-ϕ))"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,PhaseShiftGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    cosϕ2 = cos(0.5*cg.gate.ϕ[])
    sinϕ2 = sin(0.5*cg.gate.ϕ[])

    vs = (-cosϕ2*sinϕ2) .* ρv
    vs[:, 1, :] .+= (cosϕ2*sinϕ2) .* ρv[:, 4, :]
    vs[:, 2, :] .+= (sinϕ2^2)     .* ρv[:, 3, :]
    vs[:, 3, :] .-= (sinϕ2^2)     .* ρv[:, 2, :]
    vs[:, 4, :] .+= (cosϕ2*sinϕ2) .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (P(ϕ) ρ - ρ P(-ϕ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,PhaseShiftGate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2

    shift = 2^l

    sinϕ2cosϕ2 = cos(0.5*cg.gate.ϕ[]) * sin(0.5*cg.gate.ϕ[])
    sinϕ2sinϕ2 = sin(0.5*cg.gate.ϕ[]) * sin(0.5*cg.gate.ϕ[])

    @inbounds ρ.scratch .= -sinϕ2cosϕ2 .* ρ.v

    for i in 0:4^ρ.N-1
        if (i >> k) & 1 == 0
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] += sinϕ2cosϕ2 * ρ.v[i+3shift+1]
            else
                @inbounds ρ.scratch[i+1] += sinϕ2sinϕ2 * ρ.v[i+shift+1]
            end
        else
            if (i >> l) & 1 == 0
                @inbounds ρ.scratch[i+1] -= sinϕ2sinϕ2 * ρ.v[i-shift+1]
            else
                @inbounds ρ.scratch[i+1] += sinϕ2cosϕ2 * ρ.v[i-3shift+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Tailored conjugation of density matrix by SwapGate"""
function apply(ρ::DensityMatrix, cg::CircuitGate{2,SwapGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    # manual loop unrolling such that Julia's SubArray indexing works properly
    # note: not using @views function apply(...) for this method
    vs[:, 1, :, 1, :] .= ρv[:, 1, :, 1, :]
    vs[:, 1, :, 2, :] .= ρv[:, 2, :, 1, :]
    vs[:, 1, :, 3, :] .= ρv[:, 3, :, 1, :]
    vs[:, 1, :, 4, :] .= ρv[:, 4, :, 1, :]

    vs[:, 2, :, 1, :] .= ρv[:, 1, :, 2, :]
    vs[:, 2, :, 2, :] .= ρv[:, 2, :, 2, :]
    vs[:, 2, :, 3, :] .= ρv[:, 3, :, 2, :]
    vs[:, 2, :, 4, :] .= ρv[:, 4, :, 2, :]

    vs[:, 3, :, 1, :] .= ρv[:, 1, :, 3, :]
    vs[:, 3, :, 2, :] .= ρv[:, 2, :, 3, :]
    vs[:, 3, :, 3, :] .= ρv[:, 3, :, 3, :]
    vs[:, 3, :, 4, :] .= ρv[:, 4, :, 3, :]

    vs[:, 4, :, 1, :] .= ρv[:, 1, :, 4, :]
    vs[:, 4, :, 2, :] .= ρv[:, 2, :, 4, :]
    vs[:, 4, :, 3, :] .= ρv[:, 3, :, 4, :]
    vs[:, 4, :, 4, :] .= ρv[:, 4, :, 4, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by SwapGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{2,SwapGate})
    # qubit indices the gate acts on
    j, k = cg.iwire
    j, k = j < k ? (j, k) : (k, j)  # sort them
    l1 = 2j-2
    l2 = 2k-2

    shift1 = 2^l1
    shift2 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> l1) & 3
        index2 = (i >> l2) & 3

        index = i + (index1 - index2)*shift2 + (index2 - index1)*shift1
        @inbounds ρ.scratch[i+1] = ρ.v[index+1]
    end
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of 1/2 (SWAP ρ + ρ SWAP)"""
function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{2,SwapGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    vs = 0.5 .* ρv

    vs[:, 1, :, 2, :] .+= 0.5 .* ρv[:, 2, :, 1, :]
    vs[:, 1, :, 3, :] .+= 0.5 .* ρv[:, 3, :, 1, :]
    vs[:, 1, :, 4, :] .+= 0.5 .* ρv[:, 4, :, 1, :]

    vs[:, 2, :, 1, :] .+= 0.5 .* ρv[:, 1, :, 2, :]
    vs[:, 2, :, 3, :] .+= 0.5 .* ρv[:, 3, :, 2, :]
    vs[:, 2, :, 4, :] .+= 0.5 .* ρv[:, 4, :, 2, :]

    vs[:, 3, :, 1, :] .+= 0.5 .* ρv[:, 1, :, 3, :]
    vs[:, 3, :, 2, :] .+= 0.5 .* ρv[:, 2, :, 3, :]
    vs[:, 3, :, 4, :] .+= 0.5 .* ρv[:, 4, :, 3, :]

    vs[:, 4, :, 1, :] .+= 0.5 .* ρv[:, 1, :, 4, :]
    vs[:, 4, :, 2, :] .+= 0.5 .* ρv[:, 2, :, 4, :]
    vs[:, 4, :, 3, :] .+= 0.5 .* ρv[:, 3, :, 4, :]

    vs[:, 1, :, 1, :] .+= 0.5 .* (ρv[:, 2, :, 2, :] .+ ρv[:, 3, :, 3, :] .+ ρv[:, 4, :, 4, :])
    vs[:, 2, :, 2, :] .+= 0.5 .* (ρv[:, 1, :, 1, :] .- ρv[:, 3, :, 3, :] .- ρv[:, 4, :, 4, :])
    vs[:, 3, :, 3, :] .+= 0.5 .* (ρv[:, 1, :, 1, :] .- ρv[:, 2, :, 2, :] .- ρv[:, 4, :, 4, :])
    vs[:, 4, :, 4, :] .+= 0.5 .* (ρv[:, 1, :, 1, :] .- ρv[:, 2, :, 2, :] .- ρv[:, 3, :, 3, :])

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (SWAP ρ + ρ SWAP)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{2,SwapGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit indices the gate acts on
    j, k = cg.iwire
    j, k = j < k ? (j, k) : (k, j)  # sort them
    l1 = 2j-2
    l2 = 2k-2

    shift1 = 2^l1
    shift2 = 2^l2

    @inbounds ρ.scratch .= (0.5+factor) .* ρ.v

    for i in 0:4^ρ.N-1
        index1 = (i >> l1) & 3
        index2 = (i >> l2) & 3

        if index1 == index2
            if index1÷2 == 0 
                if index1%2 == 0
                    @inbounds ρ.scratch[i+1] += 0.5 * (ρ.v[i+shift1+shift2+1] + ρ.v[i+2shift1+2shift2+1] + ρ.v[i+3shift1+3shift2+1])
                else
                    @inbounds ρ.scratch[i+1] += 0.5 * (ρ.v[i-shift1-shift2+1] - ρ.v[i+shift1+shift2+1] - ρ.v[i+2shift1+2shift2+1])
                end
            else
                if index1%2 == 0
                    @inbounds ρ.scratch[i+1] += 0.5 * (ρ.v[i-2shift1-2shift2+1] - ρ.v[i-shift1-shift2+1] - ρ.v[i+shift1+shift2+1])
                else
                    @inbounds ρ.scratch[i+1] += 0.5 * (ρ.v[i-3shift1-3shift2+1] - ρ.v[i-2shift1-2shift2+1] - ρ.v[i-shift1-shift2+1])
                end
            end
        else
            index = i + (index1 - index2)*shift2 + (index2 - index1)*shift1
            @inbounds ρ.scratch[i+1] += 0.5ρ.v[index+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (SWAP ρ - ρ SWAP)"""
function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{2,SwapGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)

    vs[:, 1, :, 1, :] .= 0
    vs[:, 2, :, 2, :] .= 0
    vs[:, 3, :, 3, :] .= 0
    vs[:, 4, :, 4, :] .= 0

    vs[:, 1, :, 2, :] .= 0.5 .* (ρv[:, 4, :, 3, :] .- ρv[:, 3, :, 4, :])
    vs[:, 2, :, 1, :] .= -vs[:, 1, :, 2, :]

    vs[:, 1, :, 3, :] .= 0.5 .* (ρv[:, 2, :, 4, :] .- ρv[:, 4, :, 2, :])
    vs[:, 3, :, 1, :] .= -vs[:, 1, :, 3, :]

    vs[:, 1, :, 4, :] .= 0.5 .* (ρv[:, 3, :, 2, :] .- ρv[:, 2, :, 3, :])
    vs[:, 4, :, 1, :] .= -vs[:, 1, :, 4, :]

    vs[:, 2, :, 3, :] .= 0.5 .* (ρv[:, 1, :, 4, :] .- ρv[:, 4, :, 1, :])
    vs[:, 3, :, 2, :] .= -vs[:, 2, :, 3, :]

    vs[:, 2, :, 4, :] .= 0.5 .* (ρv[:, 3, :, 1, :] .- ρv[:, 1, :, 3, :])
    vs[:, 4, :, 2, :] .= -vs[:, 2, :, 4, :]

    vs[:, 3, :, 4, :] .= 0.5 .* (ρv[:, 1, :, 2, :] .- ρv[:, 2, :, 1, :])
    vs[:, 4, :, 3, :] .= -vs[:, 3, :, 4, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (SWAP ρ - ρ SWAP)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{2,SwapGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    l1 = 2i-2
    l2 = 2j-2

    shift1 = 2^l1
    shift2 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> l1) & 3
        index2 = (i >> l2) & 3

        if index1 == index2
            ρ.scratch[i+1] = 0
        else
            if index2 < index1
                if index2 == 0
                    if index1 == 1
                        @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i+3shift2+shift1+1] - ρ.v[i+2shift2+2shift1+1])
                    elseif index1 == 2
                        @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i+shift2+shift1+1] - ρ.v[i+3shift2-shift1+1])
                    else
                        @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i+2shift2-2shift1+1] - ρ.v[i+shift2-shift1+1])
                    end
                elseif index2 == 1
                    if index1 == 2
                        @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i-shift2+shift1+1] - ρ.v[i+2shift2-2shift1+1])                        
                    else
                        @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i+shift2-3shift1+1] - ρ.v[i-shift2-shift1+1])                        
                    end
                else
                    @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i-2shift2-2shift1+1] - ρ.v[i-shift2-3shift1+1])
                end
            else
                index = i + (index1 - index2)*shift2 + (index2 - index1)*shift1
                @inbounds ρ.scratch[i+1] = -ρ.scratch[index+1]
            end
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Tailored conjugation of density matrix by EntanglementXXGate"""
function apply(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementXXGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)

    vs[:, 1, :, 1, :] .= ρv[:, 1, :, 1, :]
    vs[:, 2, :, 2, :] .= ρv[:, 2, :, 2, :]
    vs[:, 3, :, 3, :] .= ρv[:, 3, :, 3, :]
    vs[:, 4, :, 4, :] .= ρv[:, 4, :, 4, :]

    vs[:, 1, :, 2, :] .= ρv[:, 1, :, 2, :]
    vs[:, 2, :, 1, :] .= ρv[:, 2, :, 1, :]
    vs[:, 3, :, 4, :] .= ρv[:, 3, :, 4, :]
    vs[:, 4, :, 3, :] .= ρv[:, 4, :, 3, :]

    vs[:, 1, :, 4, :] .= cosθ .* ρv[:, 1, :, 4, :] .+ sinθ .* ρv[:, 2, :, 3, :]
    vs[:, 2, :, 3, :] .= cosθ .* ρv[:, 2, :, 3, :] .- sinθ .* ρv[:, 1, :, 4, :]

    vs[:, 2, :, 4, :] .= cosθ .* ρv[:, 2, :, 4, :] .+ sinθ .* ρv[:, 1, :, 3, :]
    vs[:, 1, :, 3, :] .= cosθ .* ρv[:, 1, :, 3, :] .- sinθ .* ρv[:, 2, :, 4, :]

    vs[:, 4, :, 2, :] .= cosθ .* ρv[:, 4, :, 2, :] .+ sinθ .* ρv[:, 3, :, 1, :]
    vs[:, 3, :, 1, :] .= cosθ .* ρv[:, 3, :, 1, :] .- sinθ .* ρv[:, 4, :, 2, :]

    vs[:, 4, :, 1, :] .= cosθ .* ρv[:, 4, :, 1, :] .+ sinθ .* ρv[:, 3, :, 2, :]
    vs[:, 3, :, 2, :] .= cosθ .* ρv[:, 3, :, 2, :] .- sinθ .* ρv[:, 4, :, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by EntanglementXXGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementXXGate})

    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    l1 = 2i-2
    l2 = 2j-2

    shift1 = 2^l1
    shift2 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> l1) & 3
        index2 = (i >> l2) & 3

        if index1 == index2 || index1 ⊻ index2 == 1
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        else
            if index1 == 2 || index2 == 2
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] + sign * sinθ * ρ.v[i + (-1)^index1*shift1 + (-1)^index2*shift2+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (Rxx(θ) ρ + ρ Rxx(-θ))"""
function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementXXGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = cosθ2 .* ρv

    vs[:, 1, :, 4, :] .+= sinθ2 .* ρv[:, 2, :, 3, :]
    vs[:, 4, :, 1, :] .+= sinθ2 .* ρv[:, 3, :, 2, :]
    vs[:, 2, :, 3, :] .-= sinθ2 .* ρv[:, 1, :, 4, :]
    vs[:, 3, :, 2, :] .-= sinθ2 .* ρv[:, 4, :, 1, :]

    vs[:, 2, :, 4, :] .+= sinθ2 .* ρv[:, 1, :, 3, :]
    vs[:, 4, :, 2, :] .+= sinθ2 .* ρv[:, 3, :, 1, :]
    vs[:, 1, :, 3, :] .-= sinθ2 .* ρv[:, 2, :, 4, :]
    vs[:, 3, :, 1, :] .-= sinθ2 .* ρv[:, 4, :, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Rxx(θ) ρ + ρ Rxx(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementXXGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit indices the gate acts on

    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    l1 = 2i-2
    l2 = 2j-2

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v

    shift1 = 2^l1
    shift2 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> l1) & 3
        index2 = (i >> l2) & 3

        if index1 == index2 || index1 ⊻ index2 == 1
            continue
        else
            index = i 
            if index1%2 == 0
                index += shift1
            else
                index -= shift1
            end

            if index2%2 == 0
                index += shift2
            else
                index -= shift2
            end
            
            if index1 == 2 || index2 == 2
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] += sign * sinθ2 * ρ.v[index+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Rxx(θ) ρ - ρ Rxx(-θ))"""
function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementXXGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = zero(ρv)

    vs[:, 1, :, 1, :] .=  sinθ2 .* ρv[:, 2, :, 2, :]
    vs[:, 1, :, 2, :] .=  sinθ2 .* ρv[:, 2, :, 1, :]
    vs[:, 2, :, 1, :] .=  sinθ2 .* ρv[:, 1, :, 2, :]
    vs[:, 2, :, 2, :] .=  sinθ2 .* ρv[:, 1, :, 1, :]

    vs[:, 3, :, 3, :] .= -sinθ2 .* ρv[:, 4, :, 4, :]
    vs[:, 3, :, 4, :] .=  sinθ2 .* ρv[:, 4, :, 3, :]
    vs[:, 4, :, 3, :] .=  sinθ2 .* ρv[:, 3, :, 4, :]
    vs[:, 4, :, 4, :] .= -sinθ2 .* ρv[:, 3, :, 3, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Rxx(θ) ρ - ρ Rxx(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementXXGate})
    sinθ2 = sin(0.5*cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || index1 == index3
            index = i
            if index2 == 0
                index += shift2
            else
                index -= shift2
            end
            
            if index4 == 0
                index += shift4
            else
                index -= shift4
            end

            if a == b && index3 == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = sign * sinθ2 * ρ.v[index+1]
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Tailored conjugation of density matrix by EntanglementYYGate"""
function apply(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementYYGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)

    vs[:, 1, :, 1, :] .= ρv[:, 1, :, 1, :]
    vs[:, 2, :, 2, :] .= ρv[:, 2, :, 2, :]
    vs[:, 3, :, 3, :] .= ρv[:, 3, :, 3, :]
    vs[:, 4, :, 4, :] .= ρv[:, 4, :, 4, :]

    vs[:, 1, :, 3, :] .= ρv[:, 1, :, 3, :]
    vs[:, 3, :, 1, :] .= ρv[:, 3, :, 1, :]
    vs[:, 2, :, 4, :] .= ρv[:, 2, :, 4, :]
    vs[:, 4, :, 2, :] .= ρv[:, 4, :, 2, :]

    vs[:, 1, :, 2, :] .= cosθ .* ρv[:, 1, :, 2, :] .+ sinθ .* ρv[:, 3, :, 4, :]
    vs[:, 3, :, 4, :] .= cosθ .* ρv[:, 3, :, 4, :] .- sinθ .* ρv[:, 1, :, 2, :]

    vs[:, 2, :, 1, :] .= cosθ .* ρv[:, 2, :, 1, :] .+ sinθ .* ρv[:, 4, :, 3, :]
    vs[:, 4, :, 3, :] .= cosθ .* ρv[:, 4, :, 3, :] .- sinθ .* ρv[:, 2, :, 1, :]

    vs[:, 2, :, 3, :] .= cosθ .* ρv[:, 2, :, 3, :] .+ sinθ .* ρv[:, 4, :, 1, :]
    vs[:, 4, :, 1, :] .= cosθ .* ρv[:, 4, :, 1, :] .- sinθ .* ρv[:, 2, :, 3, :]

    vs[:, 3, :, 2, :] .= cosθ .* ρv[:, 3, :, 2, :] .+ sinθ .* ρv[:, 1, :, 4, :]
    vs[:, 1, :, 4, :] .= cosθ .* ρv[:, 1, :, 4, :] .- sinθ .* ρv[:, 3, :, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by EntanglementYYGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementYYGate})
    # qubit indices the gate acts on
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || index2 == index4
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        else
            index = i 
            if index3 == 0
                index += shift3
            else
                index -= shift3
            end

            if index1 == 0
                index += shift1
            else
                index -= shift1
            end
            
            if a == 3 || b == 3
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] + sign * sinθ * ρ.v[index+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Tailored implementation of 1/2 (Ryy(θ) ρ + ρ Ryy(-θ))"""
function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementYYGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = cosθ2 .* ρv

    vs[:, 1, :, 2, :] .+= sinθ2 .* ρv[:, 3, :, 4, :]
    vs[:, 2, :, 1, :] .+= sinθ2 .* ρv[:, 4, :, 3, :]
    vs[:, 3, :, 4, :] .-= sinθ2 .* ρv[:, 1, :, 2, :]
    vs[:, 4, :, 3, :] .-= sinθ2 .* ρv[:, 2, :, 1, :]

    vs[:, 1, :, 4, :] .-= sinθ2 .* ρv[:, 3, :, 2, :]
    vs[:, 4, :, 1, :] .-= sinθ2 .* ρv[:, 2, :, 3, :]
    vs[:, 3, :, 2, :] .+= sinθ2 .* ρv[:, 1, :, 4, :]
    vs[:, 2, :, 3, :] .+= sinθ2 .* ρv[:, 4, :, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Ryy(θ) ρ + ρ Ryy(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementYYGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit indices the gate acts on
    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || index2 == index4
            continue
        else
            index = i 
            if index3 == 0
                index += shift3
            else
                index -= shift3
            end

            if index1 == 0
                index += shift1
            else
                index -= shift1
            end
            
            if a == 3 || b == 3
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] += sign * sinθ2 * ρ.v[index+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Ryy(θ) ρ - ρ Ryy(-θ))"""
function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementYYGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = zero(ρv)

    vs[:, 1, :, 1, :] .=  sinθ2 .* ρv[:, 3, :, 3, :]
    vs[:, 1, :, 3, :] .=  sinθ2 .* ρv[:, 3, :, 1, :]
    vs[:, 3, :, 1, :] .=  sinθ2 .* ρv[:, 1, :, 3, :]
    vs[:, 3, :, 3, :] .=  sinθ2 .* ρv[:, 1, :, 1, :]

    vs[:, 2, :, 2, :] .= -sinθ2 .* ρv[:, 4, :, 4, :]
    vs[:, 2, :, 4, :] .=  sinθ2 .* ρv[:, 4, :, 2, :]
    vs[:, 4, :, 2, :] .=  sinθ2 .* ρv[:, 2, :, 4, :]
    vs[:, 4, :, 4, :] .= -sinθ2 .* ρv[:, 2, :, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Ryy(θ) ρ - ρ Ryy(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementYYGate})
    # qubit indices the gate acts on
    sinθ2 = sin(0.5*cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || index2 == index4
            index = i
            if index3 == 0
                index += shift3
            else
                index -= shift3
            end

            if index1 == 0
                index += shift1
            else
                index -= shift1
            end

            if a == b && index4 == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = sign * sinθ2 * ρ.v[index+1]
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored conjugation of density matrix by EntanglementZZGate"""
function apply(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementZZGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)

    vs[:, 1, :, 1, :] .= ρv[:, 1, :, 1, :]
    vs[:, 2, :, 2, :] .= ρv[:, 2, :, 2, :]
    vs[:, 3, :, 3, :] .= ρv[:, 3, :, 3, :]
    vs[:, 4, :, 4, :] .= ρv[:, 4, :, 4, :]

    vs[:, 1, :, 4, :] .= ρv[:, 1, :, 4, :]
    vs[:, 4, :, 1, :] .= ρv[:, 4, :, 1, :]
    vs[:, 3, :, 2, :] .= ρv[:, 3, :, 2, :]
    vs[:, 2, :, 3, :] .= ρv[:, 2, :, 3, :]

    vs[:, 1, :, 3, :] .= cosθ .* ρv[:, 1, :, 3, :] .+ sinθ .* ρv[:, 4, :, 2, :]
    vs[:, 4, :, 2, :] .= cosθ .* ρv[:, 4, :, 2, :] .- sinθ .* ρv[:, 1, :, 3, :]

    vs[:, 3, :, 1, :] .= cosθ .* ρv[:, 3, :, 1, :] .+ sinθ .* ρv[:, 2, :, 4, :]
    vs[:, 2, :, 4, :] .= cosθ .* ρv[:, 2, :, 4, :] .- sinθ .* ρv[:, 3, :, 1, :]

    vs[:, 3, :, 4, :] .= cosθ .* ρv[:, 3, :, 4, :] .+ sinθ .* ρv[:, 2, :, 1, :]
    vs[:, 2, :, 1, :] .= cosθ .* ρv[:, 2, :, 1, :] .- sinθ .* ρv[:, 3, :, 4, :]

    vs[:, 4, :, 3, :] .= cosθ .* ρv[:, 4, :, 3, :] .+ sinθ .* ρv[:, 1, :, 2, :]
    vs[:, 1, :, 2, :] .= cosθ .* ρv[:, 1, :, 2, :] .- sinθ .* ρv[:, 4, :, 3, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by EntanglementZZGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementZZGate})
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || xor(index1, index3) & xor(index2, index4) == 1
            @inbounds ρ.scratch[i+1] = ρ.v[i+1]
        else
            index = i 
            if b == 0
                index += 3shift4
            elseif b == 1
                index += shift4
            elseif b == 2
                index -= shift4
            else
                index -= 3shift4
            end

            if a == 0
                index += 3shift2
            elseif a == 1
                index += shift2
            elseif a == 2
                index -= shift2
            else
                index -= 3shift2
            end


            if a == 1 || b == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = cosθ * ρ.v[i+1] + sign * sinθ * ρ.v[index+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of 1/2 (Rzz(θ) ρ + ρ Rzz(-θ))"""
function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementZZGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = cosθ2 .* ρv

    vs[:, 1, :, 3, :] .+= sinθ2 .* ρv[:, 4, :, 2, :]
    vs[:, 3, :, 1, :] .+= sinθ2 .* ρv[:, 2, :, 4, :]
    vs[:, 4, :, 2, :] .-= sinθ2 .* ρv[:, 1, :, 3, :]
    vs[:, 2, :, 4, :] .-= sinθ2 .* ρv[:, 3, :, 1, :]

    vs[:, 4, :, 3, :] .+= sinθ2 .* ρv[:, 1, :, 2, :]
    vs[:, 3, :, 4, :] .+= sinθ2 .* ρv[:, 2, :, 1, :]
    vs[:, 1, :, 2, :] .-= sinθ2 .* ρv[:, 4, :, 3, :]
    vs[:, 2, :, 1, :] .-= sinθ2 .* ρv[:, 3, :, 4, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (Rzz(θ) ρ + ρ Rzz(-θ))"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementZZGate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit indices the gate acts on
    cosθ2 = cos(0.5*cg.gate.θ[])
    sinθ2 = sin(0.5*cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    @inbounds ρ.scratch .= (cosθ2+factor) .* ρ.v

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || xor(index1, index3) & xor(index2, index4) == 1
            continue
        else
            index = i 
            if b == 0
                index += 3shift4
            elseif b == 1
                index += shift4
            elseif b == 2
                index -= shift4
            else
                index -= 3shift4
            end

            if a == 0
                index += 3shift2
            elseif a == 1
                index += shift2
            elseif a == 2
                index -= shift2
            else
                index -= 3shift2
            end

            if a == 1 || b == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] += sign * sinθ2 * ρ.v[index+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end


"""Tailored implementation of i/2 (Rzz(θ) ρ - ρ Rzz(-θ))"""
function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementZZGate})
    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4^(i-1), 4, 4^(j-i-1), 4, 4^(ρ.N-j))

    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = zero(ρv)

    vs[:, 1, :, 1, :] .=  sinθ2 .* ρv[:, 4, :, 4, :]
    vs[:, 1, :, 4, :] .=  sinθ2 .* ρv[:, 4, :, 1, :]
    vs[:, 4, :, 1, :] .=  sinθ2 .* ρv[:, 1, :, 4, :]
    vs[:, 4, :, 4, :] .=  sinθ2 .* ρv[:, 1, :, 1, :]

    vs[:, 3, :, 3, :] .= -sinθ2 .* ρv[:, 2, :, 2, :]
    vs[:, 3, :, 2, :] .=  sinθ2 .* ρv[:, 2, :, 3, :]
    vs[:, 2, :, 3, :] .=  sinθ2 .* ρv[:, 3, :, 2, :]
    vs[:, 2, :, 2, :] .= -sinθ2 .* ρv[:, 3, :, 3, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Rzz(θ) ρ - ρ Rzz(-θ))"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{2,EntanglementZZGate})
    # qubit indices the gate acts on
    sinθ2 = sin(0.5*cg.gate.θ[])

    # qubit indices the gate acts on
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them

    k1 = 2i-1
    l1 = 2i-2
    k2 = 2j-1
    l2 = 2j-2

    shift1 = 2^k1
    shift2 = 2^l1
    shift3 = 2^k2
    shift4 = 2^l2

    for i in 0:4^ρ.N-1
        index1 = (i >> k1) & 1
        index2 = (i >> l1) & 1
        index3 = (i >> k2) & 1
        index4 = (i >> l2) & 1
        a = index1 * 2 + index2
        b = index3 * 2 + index4

        if a == b || xor(index1, index3) & xor(index2, index4) == 1
            index = i
            if b == 0
                index += 3shift4
            elseif b == 1
                index += shift4
            elseif b == 2
                index -= shift4
            else
                index -= 3shift4
            end

            if a == 0
                index += 3shift2
            elseif a == 1
                index += shift2
            elseif a == 2
                index -= shift2
            else
                index -= 3shift2
            end

            if a == b && xor(index1, index2) == 1
                sign = -1
            else
                sign = 1
            end

            @inbounds ρ.scratch[i+1] = sign * sinθ2 * ρ.v[index+1]
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end



"""Allow 'kron' to be called with a single matrix, simply returning the matrix."""
Base.kron(a::AbstractMatrix{T}) where {T} = a


"""
Projection |1><1| (auxiliary "gate" required for implementing controlled gates)
"""
struct Proj1Gate <: AbstractGate end
const proj1GateMatrix = ComplexQ[0 0; 0 1]
const sproj1GateMatrix = sparse(proj1GateMatrix)
matrix(::Proj1Gate) = proj1GateMatrix
sparse_matrix(::Proj1Gate) = sproj1GateMatrix

LinearAlgebra.ishermitian(::Proj1Gate) = true

Base.adjoint(g::Proj1Gate) = g

# wires
num_wires(::Proj1Gate)::Int = 1


"""Tailored conjugation of density matrix by Proj1Gate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{1,Proj1Gate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    # |1><1| I |1><1| =  |1><1| =  1/2 (I - Z)
    # |1><1| Z |1><1| = -|1><1| = -1/2 (I - Z)
    vs = similar(ρv)
    vs[:, 1, :] .= 0.5 .* (ρv[:, 1, :] .- ρv[:, 4, :])
    vs[:, 2, :] .= 0
    vs[:, 3, :] .= 0
    vs[:, 4, :] .= 0.5 .* (ρv[:, 4, :] .- ρv[:, 1, :])

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored conjugation of density matrix by Proj1Gate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{1,Proj1Gate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    for i in 0:4^ρ.N-1
        index1 = (i >> k) & 1
        index2 = (i >> l) & 1

        if xor(index1, index2) == 1
            @inbounds ρ.scratch[i+1] = 0
        else
            if index1 == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = 0.5 * (ρ.v[i+1] - ρ.v[i+sign*3shift+1])
        end
    end

    @inbounds ρ.v .= ρ.scratch
    return ρ
end

"""Tailored implementation of 1/2 (|1><1| ρ + ρ |1><1|)"""
@views function apply_mixed_add(ρ::DensityMatrix, cg::CircuitGate{1,Proj1Gate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = 0.5 .* ρv
    vs[:, 1, :] .-= 0.5 .* ρv[:, 4, :]
    vs[:, 4, :] .-= 0.5 .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (|1><1| ρ + ρ |1><1|)"""
@views function apply_mixed_add!(ρ::DensityMatrix, cg::CircuitGate{1,Proj1Gate}, factor::FloatQ=convert(FloatQ, 0))
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    factor = factor + 0.5
    @inbounds ρ.scratch .= factor .* ρ.v
    for i in 0:4^ρ.N-1
        index1 = (i >> k) & 1
        index2 = (i >> l) & 1

        if xor(index1, index2) == 1
            continue
        else
            if index1 == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] -= 0.5 * ρ.v[i+sign*3shift+1]
        end
    end

    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored implementation of i/2 (|1><1| ρ - ρ |1><1|)"""
@views function apply_mixed_sub(ρ::DensityMatrix, cg::CircuitGate{1,Proj1Gate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = zero(ρv)
    vs[:, 2, :] .= -0.5 .* ρv[:, 3, :]
    vs[:, 3, :] .=  0.5 .* ρv[:, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (|1><1| ρ - ρ |1><1|)"""
@views function apply_mixed_sub!(ρ::DensityMatrix, cg::CircuitGate{1,Proj1Gate})
    # qubit index the gate acts on
    j = cg.iwire[1]
    k = 2*cg.iwire[1]-1
    l = 2*cg.iwire[1]-2
    shift = 2^l

    for i in 0:4^ρ.N-1
        index1 = (i >> k) & 1
        index2 = (i >> l) & 1

        if xor(index1, index2) == 1
            if index2 == 1
                sign = -1
            else
                sign = 1
            end
            @inbounds ρ.scratch[i+1] = 0.5 * sign * ρ.v[i-sign*shift+1]
        else
            @inbounds ρ.scratch[i+1] = 0
        end
    end
    ρ.v, ρ.scratch = ρ.scratch, ρ.v
    return ρ
end

"""Tailored conjugation of density matrix by a general ControlledGate"""
function apply(ρ::DensityMatrix, cg::CircuitGate{M,ControlledGate{G}}) where {M,G}

    # number of target wires
    T = target_wires(cg.gate)

    # number of control wires
    C = control_wires(cg.gate)

    # for a single control qubit:
    # controlled-U = |0><0| x I + |1><1| x U = I + |1><1| x (U - I)
    # (generalizes to several control qubits)

    # identity term
    ρs = deepcopy(ρ)

    cgU = CircuitGate{T,G}(cg.iwire[1:T], cg.gate.U)

    # mixed term |1><1| x (U - I) ρ + ρ |1><1| x (U† - I)
    eo = binary_digits(C, 0)
    for p in 0:2^C-1
        binary_digits!(eo, p)
        τ = ρ
        for j in 1:C
            if eo[j] == 0
                τ = apply_mixed_add(τ, CircuitGate((cg.iwire[T + j],), Proj1Gate()))
            else
                τ = apply_mixed_sub(τ, CircuitGate((cg.iwire[T + j],), Proj1Gate()))
            end
        end
        nodd::Int = sum(eo)
        if iseven(nodd)
            τ.v .= apply_mixed_add(τ, cgU).v .- τ.v
        else
            τ = apply_mixed_sub(τ, cgU)
        end
        # sign factor from even powers of i
        signfac = 1 - 2*(((nodd + 1) ÷ 2) % 2)
        ρs.v .+= (2.0*signfac) .* τ.v
    end

    # conjugate by |1><1| x (U - I)
    τ = ρ
    for j in 1:C
        τ = apply(τ, CircuitGate((cg.iwire[T + j],), Proj1Gate()))
    end
    ρs.v .+= apply(τ, cgU).v .- 2.0.*apply_mixed_add(τ, cgU).v .+ τ.v

    return ρs
end

"""Tailored conjugation of density matrix by a general ControlledGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{M,ControlledGate{G}}) where {M,G}

    # number of target wires
    T = target_wires(cg.gate)

    # number of control wires
    C = control_wires(cg.gate)

    # for a single control qubit:
    # controlled-U = |0><0| x I + |1><1| x U = I + |1><1| x (U - I)
    # (generalizes to several control qubits)

    # identity term
    ρs = deepcopy(ρ.v)
    τ = deepcopy(ρ)
    cgU = CircuitGate{T,G}(cg.iwire[1:T], cg.gate.U)

    # mixed term |1><1| x (U - I) ρ + ρ |1><1| x (U† - I)
    eo = binary_digits(C, 0)
    # TODO: exploit sparsity
    for p in 0:2^C-1
        binary_digits!(eo, p)
        @inbounds τ.v .= ρs
        for j in 1:C
            if eo[j] == 0
                apply_mixed_add!(τ, CircuitGate((cg.iwire[T + j],), Proj1Gate()))
            else
                apply_mixed_sub!(τ, CircuitGate((cg.iwire[T + j],), Proj1Gate()))
            end
        end
        nodd::Int = sum(eo)
        if iseven(nodd)
            apply_mixed_add!(τ, cgU, convert(FloatQ, -1.0))
        else
            apply_mixed_sub!(τ, cgU)
        end
        # sign factor from even powers of i
        signfac = 1 - 2*(((nodd + 1) ÷ 2) % 2)
        @inbounds ρ.v .+= (2.0*signfac) .* τ.v
    end

    # conjugate by |1><1| x (U - I)
    τ.v .= ρs
    for j in 1:C
        apply!(τ, CircuitGate((cg.iwire[T + j],), Proj1Gate()))
    end
    @inbounds ρs .= τ.v
    @inbounds ρ.v .+= apply!(τ, cgU).v
    @inbounds τ.v .= ρs
    @inbounds ρ.v .-= 2.0.*apply_mixed_add!(τ, cgU, promote(-0.5)).v

    return ρ
end

"""Conjugate density matrix by a general MatrixGate"""
@views function apply(ρ::DensityMatrix, cg::CircuitGate{M,MatrixGate}) where {M}

    # Pauli matrix basis (including identity matrix)
    pauli = Matrix{ComplexQ}[Matrix{ComplexQ}(I, 2, 2), matrix(X), matrix(Y), matrix(Z)]
    halfpauli = Matrix{ComplexQ}[Matrix{ComplexQ}(0.5I, 2, 2), 0.5*matrix(X), 0.5*matrix(Y), 0.5*matrix(Z)]

    U::Matrix{ComplexQ} = matrix(cg.gate)
    # represent conjugation by U with respect to Pauli basis
    m = quaternary_digits(M, 0)
    conjU = FloatQ[real(tr(kron([pauli[p+1] for p in reverse(quaternary_digits!(m, i))]...) * U * kron([halfpauli[p+1] for p in reverse(quaternary_digits!(m, j))]...) * U'))
                for i in 0:4^M-1,
                    j in 0:4^M-1]

    ρv = reshape(ρ.v, fill(4, ρ.N)...)

    # apply conjU to circuit gate wires
    vs = similar(ρv)
    for i in 1:4^M
        # cannot use .= here since broadcasting fails for scalar numbers
        vs[sliced_index(quaternary_digits!(m, i - 1), cg.iwire, ρ.N)...] = sum(conjU[i, j] .* ρv[sliced_index(quaternary_digits!(m, j - 1), cg.iwire, ρ.N)...] for j in 1:4^M)
    end

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Conjugate density matrix by a general MatrixGate"""
@views function apply!(ρ::DensityMatrix, cg::CircuitGate{M,MatrixGate}) where {M}

    # Pauli matrix basis (including identity matrix)
    pauli = Matrix{ComplexQ}[Matrix{ComplexQ}(I, 2, 2), matrix(X), matrix(Y), matrix(Z)]
    halfpauli = Matrix{ComplexQ}[Matrix{ComplexQ}(0.5I, 2, 2), 0.5*matrix(X), 0.5*matrix(Y), 0.5*matrix(Z)]

    U::Matrix{ComplexQ} = matrix(cg.gate)
    # represent conjugation by U with respect to Pauli basis
    m = quaternary_digits(M, 0)
    conjU = FloatQ[real(tr(kron([pauli[p+1] for p in reverse(quaternary_digits!(m, i))]...) * U * kron([halfpauli[p+1] for p in reverse(quaternary_digits!(m, j))]...) * U'))
                for i in 0:4^M-1,
                    j in 0:4^M-1]

    ρv = reshape(ρ.v, fill(4, ρ.N)...)

    # apply conjU to circuit gate wires
    vs = similar(ρv)
    for i in 1:4^M
        # cannot use .= here since broadcasting fails for scalar numbers
        vs[sliced_index(quaternary_digits!(m, i - 1), cg.iwire, ρ.N)...] = sum(conjU[i, j] .* ρv[sliced_index(quaternary_digits!(m, j - 1), cg.iwire, ρ.N)...] for j in 1:4^M)
    end
    ρ.v .= reshape(vs, :)
    return ρ
end


function apply(ρ::DensityMatrix, moment::Moment) 
    for cg in moment
        ρ = apply(ρ, cg)
    end
    return ρ
end

function apply!(ρ::DensityMatrix, moment::Moment) 
    for cg in moment
        apply!(ρ, cg)
    end
    return ρ
end

function apply(ρ::DensityMatrix, moments::Vector{Moment})
    for m in moments
        ρ = apply(ρ, m)
    end
    return ρ
end

function apply!(ρ::DensityMatrix, moments::Vector{Moment})
    for m in moments
        apply!(ρ, m)
    end
    return ρ
end

"""
    apply(ρ::DensityMatrix, c::Circuit{N}) where {N}

Compute list of expectation values from measurement operators in `c.meas`, after applying circuit gates in `c.cgc` on the N-qubit density matrix `ρ`
"""
function apply(ρ::DensityMatrix, c::Circuit{N}) where {N}
    ρ.N == N || error("Qubit number $N of circuit not equal to qubit number $(ρ.N) of density matrix")
    ρ = apply(ρ, c.moments)
    ρmat = matrix(ρ)
    return [real(tr(sparse_matrix(m, N) * ρmat)) for m in c.meas]
end


"""
    apply(ρ::DensityMatrix, c::Circuit{N}) where {N}

Compute list of expectation values from measurement operators in `c.meas`, after applying circuit gates in `c.cgc` on the N-qubit density matrix `ρ`
"""
function apply!(ρ::DensityMatrix, c::Circuit{N}) where {N}
    ρ.N == N || error("Qubit number $N of circuit not equal to qubit number $(ρ.N) of density matrix")
    apply!(ρ, c.moments)
    ρmat = matrix(ρ)
    return [real(tr(sparse_matrix(m, N) * ρmat)) for m in c.meas]
end

"""
    apply(ρ::DensityMatrix, cgs::Vector{<:CircuitGate})

Apply a [`CircuitGate`](@ref) to a quantum density matrix `ρ`.
"""
function apply(ρ::DensityMatrix, cgs::Vector{<:CircuitGate})
    req = maximum(req_wires.(cgs))
    ρ.N >= req || error("Qubit number $req of [CircuitGate](@ref) not equal to qubit number $(ρ.N) of density matrix `ρ`")
    for cg in cgs
        ρ = apply(ρ, cg)
    end
    return ρ
end

"""
    apply(ρ::DensityMatrix, cgs::Vector{<:CircuitGate})

Apply a [`CircuitGate`](@ref) to a quantum density matrix `ρ`.
"""
function apply!(ρ::DensityMatrix, cgs::Vector{<:CircuitGate})
    req = maximum(req_wires.(cgs))
    ρ.N >= req || error("Qubit number $req of [CircuitGate](@ref) not equal to qubit number $(ρ.N) of density matrix `ρ`")
    for cg in cgs
        apply!(ρ, cg)
    end
    return ρ
end