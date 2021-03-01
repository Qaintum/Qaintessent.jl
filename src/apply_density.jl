

"""Tailored conjugation of density matrix by XGate"""
@views function apply(cg::CircuitGate{1,XGate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # X I X =  I
    vs[:, 2, :] .=  ρv[:, 2, :]     # X X X =  X
    vs[:, 3, :] .= -ρv[:, 3, :]     # X Y X = -Y
    vs[:, 4, :] .= -ρv[:, 4, :]     # X Z X = -Z

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of 1/2 (X ρ + ρ X)"""
@views function apply_mixed_add(cg::CircuitGate{1,XGate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = similar(ρv)
    vs[:, 1, :] .= ρv[:, 2, :]      # 1/2 (X X + X X) = I
    vs[:, 2, :] .= ρv[:, 1, :]      # 1/2 (X I + I X) = X
    vs[:, 3, :] .= 0                # 1/2 (X Y + Y X) = 0
    vs[:, 4, :] .= 0                # 1/2 (X Z + Z X) = 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (X ρ - ρ X)"""
@views function apply_mixed_sub(cg::CircuitGate{1,XGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by YGate"""
@views function apply(cg::CircuitGate{1,YGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Y ρ + ρ Y)"""
@views function apply_mixed_add(cg::CircuitGate{1,YGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Y ρ - ρ Y)"""
@views function apply_mixed_sub(cg::CircuitGate{1,YGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by ZGate"""
@views function apply(cg::CircuitGate{1,ZGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Z ρ + ρ Z)"""
@views function apply_mixed_add(cg::CircuitGate{1,ZGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Z ρ - ρ Z)"""
@views function apply_mixed_sub(cg::CircuitGate{1,ZGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by the Hadamard gate"""
@views function apply(cg::CircuitGate{1,HadamardGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (H ρ + ρ H)"""
@views function apply_mixed_add(cg::CircuitGate{1,HadamardGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (H ρ - ρ H)"""
@views function apply_mixed_sub(cg::CircuitGate{1,HadamardGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by SGate"""
@views function apply(cg::CircuitGate{1,SGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (S ρ + ρ S†)"""
@views function apply_mixed_add(cg::CircuitGate{1,SGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (S ρ - ρ S†)"""
@views function apply_mixed_sub(cg::CircuitGate{1,SGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by SdagGate"""
@views function apply(cg::CircuitGate{1,SdagGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (S† ρ + ρ S)"""
@views function apply_mixed_add(cg::CircuitGate{1,SdagGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (S† ρ - ρ S)"""
@views function apply_mixed_sub(cg::CircuitGate{1,SdagGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by TGate"""
@views function apply(cg::CircuitGate{1,TGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (T ρ + ρ T†)"""
@views function apply_mixed_add(cg::CircuitGate{1,TGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (T ρ - ρ T†)"""
@views function apply_mixed_sub(cg::CircuitGate{1,TGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by TdagGate"""
@views function apply(cg::CircuitGate{1,TdagGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (T† ρ + ρ T)"""
@views function apply_mixed_add(cg::CircuitGate{1,TdagGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (T† ρ - ρ T)"""
@views function apply_mixed_sub(cg::CircuitGate{1,TdagGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by RxGate"""
@views function apply(cg::CircuitGate{1,RxGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Rx(θ) ρ + ρ Rx(-θ))"""
@views function apply_mixed_add(cg::CircuitGate{1,RxGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Rx(θ) ρ - ρ Rx(-θ))"""
@views function apply_mixed_sub(cg::CircuitGate{1,RxGate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = similar(ρv)
    vs[:, [1, 2], :] .= sinθ2 .* ρv[:, [2, 1], :]
    vs[:, [3, 4], :] .= 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end


"""Tailored conjugation of density matrix by RyGate"""
@views function apply(cg::CircuitGate{1,RyGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Ry(θ) ρ + ρ Ry(-θ))"""
@views function apply_mixed_add(cg::CircuitGate{1,RyGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Ry(θ) ρ - ρ Ry(-θ))"""
@views function apply_mixed_sub(cg::CircuitGate{1,RyGate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = similar(ρv)
    vs[:, [1, 3], :] .= sinθ2 .* ρv[:, [3, 1], :]
    vs[:, [2, 4], :] .= 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end


"""Tailored conjugation of density matrix by RzGate"""
@views function apply(cg::CircuitGate{1,RzGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Rz(θ) ρ + ρ Rz(-θ))"""
@views function apply_mixed_add(cg::CircuitGate{1,RzGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Rz(θ) ρ - ρ Rz(-θ))"""
@views function apply_mixed_sub(cg::CircuitGate{1,RzGate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))
    sinθ2 = sin(0.5*cg.gate.θ[])

    vs = similar(ρv)
    vs[:, [1, 4], :] .= sinθ2 .* ρv[:, [4, 1], :]
    vs[:, [2, 3], :] .= 0

    return DensityMatrix(reshape(vs, :), ρ.N)
end


"""Tailored conjugation of density matrix by RotationGate"""
@views function apply(cg::CircuitGate{1,RotationGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Rn(θ) ρ + ρ Rn(-θ))"""
@views function apply_mixed_add(cg::CircuitGate{1,RotationGate}, ρ::DensityMatrix)
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
    sn = sinθ2 * n

    vs = cosθ2 .* ρv
    vs[:, 2, :] .+= (sn[2] .* ρv[:, 4, :] .- sn[3] .* ρv[:, 3, :])
    vs[:, 3, :] .+= (sn[3] .* ρv[:, 2, :] .- sn[1] .* ρv[:, 4, :])
    vs[:, 4, :] .+= (sn[1] .* ρv[:, 3, :] .- sn[2] .* ρv[:, 2, :])

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (Rn(θ) ρ - ρ Rn(-θ))"""
@views function apply_mixed_sub(cg::CircuitGate{1,RotationGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by PhaseShiftGate"""
@views function apply(cg::CircuitGate{1,PhaseShiftGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (P(ϕ) ρ + ρ P(-ϕ))"""
@views function apply_mixed_add(cg::CircuitGate{1,PhaseShiftGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (P(ϕ) ρ - ρ P(-ϕ))"""
@views function apply_mixed_sub(cg::CircuitGate{1,PhaseShiftGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by SwapGate"""
function apply(cg::CircuitGate{2,SwapGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (SWAP ρ + ρ SWAP)"""
function apply_mixed_add(cg::CircuitGate{2,SwapGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (SWAP ρ - ρ SWAP)"""
function apply_mixed_sub(cg::CircuitGate{2,SwapGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by EntanglementXXGate"""
function apply(cg::CircuitGate{2,EntanglementXXGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Rxx(θ) ρ + ρ Rxx(-θ))"""
function apply_mixed_add(cg::CircuitGate{2,EntanglementXXGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Rxx(θ) ρ - ρ Rxx(-θ))"""
function apply_mixed_sub(cg::CircuitGate{2,EntanglementXXGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by EntanglementYYGate"""
function apply(cg::CircuitGate{2,EntanglementYYGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Ryy(θ) ρ + ρ Ryy(-θ))"""
function apply_mixed_add(cg::CircuitGate{2,EntanglementYYGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Ryy(θ) ρ - ρ Ryy(-θ))"""
function apply_mixed_sub(cg::CircuitGate{2,EntanglementYYGate}, ρ::DensityMatrix)
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


"""Tailored conjugation of density matrix by EntanglementZZGate"""
function apply(cg::CircuitGate{2,EntanglementZZGate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (Rzz(θ) ρ + ρ Rzz(-θ))"""
function apply_mixed_add(cg::CircuitGate{2,EntanglementZZGate}, ρ::DensityMatrix)
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

"""Tailored implementation of i/2 (Rzz(θ) ρ - ρ Rzz(-θ))"""
function apply_mixed_sub(cg::CircuitGate{2,EntanglementZZGate}, ρ::DensityMatrix)
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


"""Allow 'kron' to be called with a single matrix, simply returning the matrix."""
Base.kron(a::AbstractMatrix{T}) where {T} = a


"""
Projection |1><1| (auxiliary "gate" required for implementing controlled gates)
"""
struct Proj1Gate <: AbstractGate end

matrix(::Proj1Gate) = ComplexF64[0 0; 0 1]
sparse_matrix(::Proj1Gate) = sparse(matrix(Proj1Gate()))

LinearAlgebra.ishermitian(::Proj1Gate) = true

Base.adjoint(g::Proj1Gate) = g

# wires
num_wires(::Proj1Gate)::Int = 1


"""Tailored conjugation of density matrix by Proj1Gate"""
@views function apply(cg::CircuitGate{1,Proj1Gate}, ρ::DensityMatrix)
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

"""Tailored implementation of 1/2 (|1><1| ρ + ρ |1><1|)"""
@views function apply_mixed_add(cg::CircuitGate{1,Proj1Gate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = 0.5 .* ρv
    vs[:, 1, :] .-= 0.5 .* ρv[:, 4, :]
    vs[:, 4, :] .-= 0.5 .* ρv[:, 1, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end

"""Tailored implementation of i/2 (|1><1| ρ - ρ |1><1|)"""
@views function apply_mixed_sub(cg::CircuitGate{1,Proj1Gate}, ρ::DensityMatrix)
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4^(j-1), 4, 4^(ρ.N-j))

    vs = zero(ρv)
    vs[:, 2, :] .= -0.5 .* ρv[:, 3, :]
    vs[:, 3, :] .=  0.5 .* ρv[:, 2, :]

    return DensityMatrix(reshape(vs, :), ρ.N)
end


"""Tailored conjugation of density matrix by a general ControlledGate"""
function apply(cg::CircuitGate{M,ControlledGate{G}}, ρ::DensityMatrix) where {M,G}

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
    for p in 0:2^C-1
        eo = binary_digits(p, C)
        τ = ρ
        for j in 1:C
            if eo[j] == 0
                τ = apply_mixed_add(CircuitGate((cg.iwire[T + j],), Proj1Gate()), τ)
            else
                τ = apply_mixed_sub(CircuitGate((cg.iwire[T + j],), Proj1Gate()), τ)
            end
        end
        nodd::Int = sum(eo)
        if iseven(nodd)
            τ.v .= apply_mixed_add(cgU, τ).v .- τ.v
        else
            τ = apply_mixed_sub(cgU, τ)
        end
        # sign factor from even powers of i
        signfac = 1 - 2*(((nodd + 1) ÷ 2) % 2)
        ρs.v .+= (2.0*signfac) .* τ.v
    end

    # conjugate by |1><1| x (U - I)
    τ = ρ
    for j in 1:C
        τ = apply(CircuitGate((cg.iwire[T + j],), Proj1Gate()), τ)
    end
    ρs.v .+= apply(cgU, τ).v .- 2.0.*apply_mixed_add(cgU, τ).v .+ τ.v

    return ρs
end


"""Conjugate density matrix by a general MatrixGate"""
@views function apply(cg::CircuitGate{M,MatrixGate}, ρ::DensityMatrix) where {M}

    # Pauli matrix basis (including identity matrix)
    pauli = Matrix{ComplexF64}[Matrix{ComplexF64}(I, 2, 2), matrix(X), matrix(Y), matrix(Z)]
    halfpauli = Matrix{ComplexF64}[Matrix{ComplexF64}(0.5I, 2, 2), 0.5*matrix(X), 0.5*matrix(Y), 0.5*matrix(Z)]

    U::Matrix{ComplexF64} = matrix(cg.gate)
    # represent conjugation by U with respect to Pauli basis
    conjU = Float64[real(tr(kron([pauli[p+1] for p in reverse(quaternary_digits(i, M))]...) * U * kron([halfpauli[p+1] for p in reverse(quaternary_digits(j, M))]...) * U'))
                for i in 0:4^M-1,
                    j in 0:4^M-1]

    ρv = reshape(ρ.v, fill(4, ρ.N)...)

    # apply conjU to circuit gate wires
    vs = similar(ρv)
    for i in 1:4^M
        # cannot use .= here since broadcasting fails for scalar numbers
        vs[sliced_index(quaternary_digits(i - 1, M), cg.iwire, ρ.N)...] = sum(conjU[i, j] .* ρv[sliced_index(quaternary_digits(j - 1, M), cg.iwire, ρ.N)...] for j in 1:4^M)
    end

    return DensityMatrix(reshape(vs, :), ρ.N)
end


function apply(moment::Moment, ρ::DensityMatrix) 
    for cg in moment
        ρ = apply(cg, ρ)
    end
    return ρ
end

function apply(moments::Vector{Moment}, ρ::DensityMatrix)
    for m in moments
        ρ = apply(m, ρ)
    end
    return ρ
end


"""
    apply(c::Circuit{N}, ρ::DensityMatrix) where {N}

Compute list of expectation values from measurement operators in `c.meas`, after applying circuit gates in `c.cgc` on the N-qubit density matrix `ρ`
"""
function apply(c::Circuit{N}, ρ::DensityMatrix) where {N}
    ρ.N == N || error("Qubit number $N of circuit not equal to qubit number $(ρ.N) of density matrix")
    ρ = apply(c.moments, ρ)
    ρmat = matrix(ρ)
    return [real(tr(sparse_matrix(m, N) * ρmat)) for m in c.meas]
end

"""
    apply(cgs::Vector{<:CircuitGate}, ρ::DensityMatrix)

Apply a [`CircuitGate`](@ref) to a quantum density matrix `ρ`.
"""
function apply(cgs::Vector{<:CircuitGate}, ρ::DensityMatrix)
    req = maximum(req_wires.(cgs))
    ρ.N >= req || error("Qubit number $req of [CircuitGate](@ref) not equal to qubit number $(ρ.N) of density matrix `ρ`")
    for cg in cgs
        ρ = apply(cg, ρ)
    end
    return ρ
end