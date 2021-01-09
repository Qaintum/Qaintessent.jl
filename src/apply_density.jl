

"""Tailored apply to density matrix for XGate"""
function apply(cg::CircuitGate{1,XGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # X I X =  I
    vs[:, 2, :] .=  ρv[:, 2, :]     # X X X =  X
    vs[:, 3, :] .= -ρv[:, 3, :]     # X Y X = -Y
    vs[:, 4, :] .= -ρv[:, 4, :]     # X Z X = -Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for YGate"""
function apply(cg::CircuitGate{1,YGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # Y I Y =  I
    vs[:, 2, :] .= -ρv[:, 2, :]     # Y X Y = -X
    vs[:, 3, :] .=  ρv[:, 3, :]     # Y Y Y =  Y
    vs[:, 4, :] .= -ρv[:, 4, :]     # Y Z Y = -Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for ZGate"""
function apply(cg::CircuitGate{1,ZGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # Z I Z =  I
    vs[:, 2, :] .= -ρv[:, 2, :]     # Z X Z = -X
    vs[:, 3, :] .= -ρv[:, 3, :]     # Z Y Z = -Y
    vs[:, 4, :] .=  ρv[:, 4, :]     # Z Z Z =  Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for HadamardGate"""
function apply(cg::CircuitGate{1,HadamardGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # H I H =  I
    vs[:, 2, :] .=  ρv[:, 4, :]     # H Z H =  X
    vs[:, 3, :] .= -ρv[:, 3, :]     # H Y H = -Y
    vs[:, 4, :] .=  ρv[:, 2, :]     # H X H =  Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for SGate"""
function apply(cg::CircuitGate{1,SGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # S I S^† =  I
    vs[:, 2, :] .= -ρv[:, 3, :]     # S Y S^† = -X
    vs[:, 3, :] .=  ρv[:, 2, :]     # S X S^† =  Y
    vs[:, 4, :] .=  ρv[:, 4, :]     # S Z S^† =  Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for SdagGate"""
function apply(cg::CircuitGate{1,SdagGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]     # S^† I S =  I
    vs[:, 2, :] .=  ρv[:, 3, :]     # S^† Y S =  X
    vs[:, 3, :] .= -ρv[:, 2, :]     # S^† X S = -Y
    vs[:, 4, :] .=  ρv[:, 4, :]     # S^† Z S =  Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for TGate"""
function apply(cg::CircuitGate{1,TGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=  ρv[:, 1, :]                         # T I T^† = I
    vs[:, 2, :] .= (ρv[:, 2, :] .- ρv[:, 3, :])/sqrt(2) # (T X T^† - T Y T^†)/sqrt(2) = X
    vs[:, 3, :] .= (ρv[:, 2, :] .+ ρv[:, 3, :])/sqrt(2) # (T X T^† + T Y T^†)/sqrt(2) = Y
    vs[:, 4, :] .=  ρv[:, 4, :]                         # T Z T^† = Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for TdagGate"""
function apply(cg::CircuitGate{1,TdagGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .=   ρv[:, 1, :]                            # T^† I T = I
    vs[:, 2, :] .= ( ρv[:, 2, :] .+ ρv[:, 3, :])/sqrt(2)    # ( T^† X T + T^† Y T)/sqrt(2) = X
    vs[:, 3, :] .= (-ρv[:, 2, :] .+ ρv[:, 3, :])/sqrt(2)    # (-T^† X T + T^† Y T)/sqrt(2) = Y
    vs[:, 4, :] .=   ρv[:, 4, :]                            # T^† Z T = Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for RxGate"""
function apply(cg::CircuitGate{1,RxGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # Rx(θ) I Rx(-θ) = I
    vs[:, 2, :] .=       ρv[:, 2, :]                        # Rx(θ) X Rx(-θ) = X
    vs[:, 3, :] .= cosθ.*ρv[:, 3, :] .- sinθ.*ρv[:, 4, :]   # cos(θ) Rx(θ) Y Rx(-θ) - sin(θ) Rx(θ) Z Rx(-θ) = Y
    vs[:, 4, :] .= sinθ.*ρv[:, 3, :] .+ cosθ.*ρv[:, 4, :]   # sin(θ) Rx(θ) Y Rx(-θ) + cos(θ) Rx(θ) Z Rx(-θ) = Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for RyGate"""
function apply(cg::CircuitGate{1,RyGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # Ry(θ) I Ry(-θ) = I
    vs[:, 2, :] .= sinθ.*ρv[:, 4, :] .+ cosθ.*ρv[:, 2, :]   # sin(θ) Ry(θ) Z Ry(-θ) + cos(θ) Ry(θ) X Ry(-θ) = X
    vs[:, 3, :] .=       ρv[:, 3, :]                        # Ry(θ) Y Ry(-θ) = Y
    vs[:, 4, :] .= cosθ.*ρv[:, 4, :] .- sinθ.*ρv[:, 2, :]   # cos(θ) Ry(θ) Z Ry(-θ) - sin(θ) Ry(θ) X Ry(-θ) = Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for RzGate"""
function apply(cg::CircuitGate{1,RzGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)
    cosθ = cos(cg.gate.θ[])
    sinθ = sin(cg.gate.θ[])

    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # Rz(θ) I Rz(-θ) = I
    vs[:, 2, :] .= cosθ.*ρv[:, 2, :] .- sinθ.*ρv[:, 3, :]   # cos(θ) Rz(θ) X Rz(-θ) - sin(θ) Rz(θ) Y Rz(-θ) = X
    vs[:, 3, :] .= sinθ.*ρv[:, 2, :] .+ cosθ.*ρv[:, 3, :]   # sin(θ) Rz(θ) X Rz(-θ) + cos(θ) Rz(θ) Y Rz(-θ) = Y
    vs[:, 4, :] .=       ρv[:, 4, :]                        # Rz(θ) Z Rz(-θ) = Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for RotationGate"""
function apply(cg::CircuitGate{1,RotationGate}, ρ::DensityMatrix{N}) where {N}
    θ = norm(cg.gate.nθ)
    if θ == 0
        # for consistency, we return a copy here
        return DensityMatrix{N}(copy(ρ.v))
    end
    cosθ = cos(θ)
    sinθ = sin(θ)
    n = cg.gate.nθ/θ

    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)

    vs = similar(ρv)
    vs[:, 1, :] .= ρv[:, 1, :]      # Rn(θ) I Rn(-θ) = I
    # Rodrigues' rotation formula
    vs[:, 2, :] .= cosθ.*ρv[:, 2, :] .+ sinθ.*(n[2].*ρv[:, 4, :] .- n[3].*ρv[:, 3, :]) .+ ((1 - cosθ)*n[1]).*(n[1].*ρv[:, 2, :] .+ n[2].*ρv[:, 3, :] .+ n[3].*ρv[:, 4, :])
    vs[:, 3, :] .= cosθ.*ρv[:, 3, :] .+ sinθ.*(n[3].*ρv[:, 2, :] .- n[1].*ρv[:, 4, :]) .+ ((1 - cosθ)*n[2]).*(n[1].*ρv[:, 2, :] .+ n[2].*ρv[:, 3, :] .+ n[3].*ρv[:, 4, :])
    vs[:, 4, :] .= cosθ.*ρv[:, 4, :] .+ sinθ.*(n[1].*ρv[:, 3, :] .- n[2].*ρv[:, 2, :]) .+ ((1 - cosθ)*n[3]).*(n[1].*ρv[:, 2, :] .+ n[2].*ρv[:, 3, :] .+ n[3].*ρv[:, 4, :])

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for PhaseShiftGate"""
function apply(cg::CircuitGate{1,PhaseShiftGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    ρv = reshape(ρ.v, 4 .^[j-1, 1, N-j]...)
    cosϕ = cos(cg.gate.ϕ[])
    sinϕ = sin(cg.gate.ϕ[])

    # agrees with Rz(θ) since global prefactor cancels
    vs = similar(ρv)
    vs[:, 1, :] .=       ρv[:, 1, :]                        # P(ϕ) I P(-ϕ) = I
    vs[:, 2, :] .= cosϕ.*ρv[:, 2, :] .- sinϕ.*ρv[:, 3, :]   # cos(ϕ) P(ϕ) X P(-ϕ) - sin(ϕ) P(ϕ) Y P(-ϕ) = X
    vs[:, 3, :] .= sinϕ.*ρv[:, 2, :] .+ cosϕ.*ρv[:, 3, :]   # sin(ϕ) P(ϕ) X P(-ϕ) + cos(ϕ) P(ϕ) Y P(-ϕ) = Y
    vs[:, 4, :] .=       ρv[:, 4, :]                        # P(ϕ) Z P(-ϕ) = Z

    return DensityMatrix{N}(reshape(vs, :))
end


"""Tailored apply to density matrix for SwapGate"""
function apply(cg::CircuitGate{2,SwapGate}, ρ::DensityMatrix{N}) where {N}

    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i)  # sort them
    ρv = reshape(ρ.v, 4 .^[i-1, 1, j-i-1, 1, N-j]...)

    vs = similar(ρv)
    for i in 1:4
        for j in 1:4
            vs[:, i, :, j, :] .= ρv[:, j, :, i, :]
        end
    end

    return DensityMatrix{N}(reshape(vs, :))
end


"""Allow 'kron' to be called with a single matrix, simply returning the matrix."""
Base.kron(a::AbstractMatrix{T}) where {T} = a


"""Tailored apply to density matrix for a general ControlledGate"""
function apply(cg::CircuitGate{M,ControlledGate{G}}, ρ::DensityMatrix{N}) where {G,N,M}

    # Pauli matrix basis (including identity matrix)
    pauli = [Matrix{Float64}(I, 2, 2), matrix(X), matrix(Y), matrix(Z)]
    halfpauli = [Matrix{Float64}(0.5I, 2, 2), 0.5*matrix(X), 0.5*matrix(Y), 0.5*matrix(Z)]
    
    # M is the length of iwire, T is the number of target wires
    T = target_wires(cg.gate)

    # number of control wires
    C = control_wires(cg.gate)

    ttuples = reshape(Qaintessent.cartesian_tuples(4, T), :)
    ctuples = reshape(Qaintessent.cartesian_tuples(4, C), :)

    # conjugation by |1><1| represented with respect to Pauli basis
    conj1X1 = [
         0.5   0     0    -0.5;
         0     0     0     0  ;
         0     0     0     0  ;
        -0.5   0     0     0.5]

    # tr(σi |1><1| σj I) for all i, j
    mult_1X1_I = [
         0.5   0     0    -0.5;
         0     0.5   0.5im 0  ;
         0    -0.5im 0.5   0  ;
        -0.5   0     0     0.5]

    U = matrix(cg.gate.U)

    # represent conjugation by (U - I) with respect to Pauli basis
    conjUI = [real(tr(kron([pauli[p+1] for p in reverse(it)]...) * (U - I) * kron([halfpauli[p+1] for p in reverse(jt)]...) * (U' - I)))
                for it in ttuples,
                    jt in ttuples]

    # represent (U - I) with respect to Pauli basis
    UI = [tr(kron([halfpauli[p+1] for p in reverse(it)]...) * (U - I)) for it in ttuples]

    # pairwise Pauli matrix multiplication phase factor table
    pauli_mult_phase = [
        1   1   1   1 ;
        1   1   im -im;
        1  -im  1   im;
        1   im -im  1 ]

    # represent multiplication by (U - I) from the left and I from the right with respect to Pauli basis;
    # bitwise XOR gives index of Pauli matrix resulting from product of two Pauli matrices
    mult_UI_I = [prod(pauli_mult_phase[jt[k]+1, it[k]+1] for k in 1:T) * UI[sum((jt[k] ⊻ it[k]) << 2(k-1) for k in 1:T) + 1]
                for it in ttuples,
                    jt in ttuples]

    ρv = reshape(ρ.v, fill(4, N)...)

    # for a single control qubit:
    # controlled-U = |0><0| x I + |1><1| x U = I + |1><1| x (U - I)
    # (generalizes to several control qubits)
    # conjugation by controlled-U represented with respect to Pauli basis:
    # I + 2*real(kron(kron(fill(mult_1X1_I, C)...), mult_UI_I)) + kron(kron(fill(conj1X1, C)...), conjUI)

    vs = copy(ρv)   # realize identity map in decomposition of controlled-U
    for (ic, icu) in enumerate(ctuples), (it, itu) in enumerate(ttuples)
        for (jc, jcu) in enumerate(ctuples), (jt, jtu) in enumerate(ttuples)
            # 'prod' realizes Kronecker product
            conj_cU = 2*real(prod(mult_1X1_I[icu[k]+1, jcu[k]+1] for k in 1:C) * mult_UI_I[it, jt]) + prod(conj1X1[icu[k]+1, jcu[k]+1] for k in 1:C) * conjUI[it, jt]
            if conj_cU != 0
                # cannot use .= here since broadcasting fails for scalar numbers
                vs[sliced_index((itu..., icu...), cg.iwire, N)...] += conj_cU .* ρv[sliced_index((jtu..., jcu...), cg.iwire, N)...]
            end
        end
    end

    return DensityMatrix{N}(reshape(vs, :))
end


"""Apply a general multiple qubit gate to a density matrix (in particular covering MatrixGate)"""
function apply(cg::CircuitGate{M,<:AbstractGate}, ρ::DensityMatrix{N}) where {M,N}

    # Pauli matrix basis (including identity matrix)
    pauli = [Matrix{Float64}(I, 2, 2), matrix(X), matrix(Y), matrix(Z)]
    halfpauli = [Matrix{Float64}(0.5I, 2, 2), 0.5*matrix(X), 0.5*matrix(Y), 0.5*matrix(Z)]

    gtuples = reshape(Qaintessent.cartesian_tuples(4, M), :)

    U = matrix(cg.gate)
    # represent conjugation by U with respect to Pauli basis
    conjU = [real(tr(kron([pauli[p+1] for p in reverse(it)]...) * U * kron([halfpauli[p+1] for p in reverse(jt)]...) * U'))
                for it in gtuples,
                    jt in gtuples]

    ρv = reshape(ρ.v, fill(4, N)...)

    # apply conjU to circuit gate wires
    vs = similar(ρv)
    for (i, it) in enumerate(gtuples)
        # cannot use .= here since broadcasting fails for scalar numbers
        vs[sliced_index(it, cg.iwire, N)...] = sum(conjU[i, j] .* ρv[sliced_index(jt, cg.iwire, N)...] for (j, jt) in enumerate(gtuples))
    end

    return DensityMatrix{N}(reshape(vs, :))
end
