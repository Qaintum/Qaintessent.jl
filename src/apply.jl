
"""Construct sliced index via target map (utility function)"""
function sliced_index(idx::Tuple, targetwires::Tuple, N::Int)
    islice = Vector{Any}(fill(Colon(), N))
    M = length(targetwires)
    for k in 1:M
        islice[targetwires[k]] = idx[k] + 1
    end
    return islice
end

"""
    apply(m::AbstractMatrix, ψ::AbstractVector)

"""
function apply(m::AbstractMatrix, ψ::Vector{<:Complex})
    return m*ψ
end

function _apply(m::AbstractMatrix, ψ::Vector{<:Complex})
    return m*ψ
end

"""
    apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}

Apply a `CircuitGate{M,G}` to a quantum state vector `ψ`.
"""
function apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")

    _apply(cg, ψ, N)
end

"""
    apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})

Apply a `CircuitGate{M,G}` to a quantum state vector `ψ`.
"""
function apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    ψl = deepcopy(ψ)
    l == 2^N || error("Vector length must be a power of 2")
    for cg in cgs 
        ψl = _apply(cg, ψl, N)
    end
    ψl
end

function _apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}, N::Int) where {M,G}
    gtuples = reshape(Qaintessent.cartesian_tuples(2, M), :)
    U = matrix(cg.gate)
    ψ = reshape(ψ, fill(2, N)...)

    ψs = similar(ψ)
    for (i, it) in enumerate(gtuples)
        # cannot use .= here since broadcasting fails for scalar numbers
        ψs[sliced_index(it, cg.iwire, N)...] = sum(U[i, j] .* ψ[sliced_index(jt, cg.iwire, N)...] for (j, jt) in enumerate(gtuples))
    end
    return ComplexF64[ψs...]
end



"""Tailored apply for XGate"""
function _apply(cg::CircuitGate{1,XGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))

    return reshape(ψ[:, [2, 1], :], :)
end


"""Tailored apply for YGate"""
function _apply(cg::CircuitGate{1,YGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)

    A[:, 1, :] .= -im.*ψ[:, 2, :]
    A[:, 2, :] .=  im.*ψ[:, 1, :]
    return reshape(A, :)
end


"""Tailored apply for ZGate"""
function _apply(cg::CircuitGate{1,ZGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .=   ψ[:, 1, :]
    A[:, 2, :] .= .-ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for HadamardGate"""
function _apply(cg::CircuitGate{1,HadamardGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= (ψ[:, 1, :] .+ ψ[:, 2, :]) ./ sqrt(2)
    A[:, 2, :] .= (ψ[:, 1, :] .- ψ[:, 2, :]) ./ sqrt(2)
    return reshape(A, :)
end


"""Tailored apply for SGate"""
function _apply(cg::CircuitGate{1,SGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= im .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for SdagGate"""
function _apply(cg::CircuitGate{1,SdagGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= -im .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for TGate"""
function _apply(cg::CircuitGate{1,TGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= Base.exp(im*π/4) .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for TdagGate"""
function _apply(cg::CircuitGate{1,TdagGate}, ψ::Vector{<:Complex}, N::Int)
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= Base.exp(-im*π/4) .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for RzGate"""
function _apply(cg::CircuitGate{1,RzGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    θ = cg.gate.θ[]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= Base.exp(-im*θ/2).*ψ[:, 1, :]
    A[:, 2, :] .= Base.exp( im*θ/2).*ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for PhaseShiftGate"""
function _apply(cg::CircuitGate{1,PhaseShiftGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= Base.exp(im*cg.gate.ϕ[]) .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for a general single qubit gate"""
function _apply(cg::CircuitGate{1,AbstractGate}, ψ::Vector{<:Complex}, N::Int)
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    U = matrix(cg.gate)
    A[:, 1, :] .= U[1, 1] .* ψ[:, 1, :] .+ U[1, 2] .* ψ[:, 2, :]
    A[:, 2, :] .= U[2, 1] .* ψ[:, 1, :] .+ U[2, 2] .* ψ[:, 2, :]
    return reshape(A, :)
end
#

"""Tailored apply for SwapGate"""
function _apply(cg::CircuitGate{2,SwapGate}, ψ::Vector{<:Complex}, N::Int)
    i,j = cg.iwire
    i,j = i<j ? (i,j) : (j,i) #sort them
    qubit_slices = [i-1, 1, j-i-1, 1, N-j]

    ψr = reshape(ψ, 2 .^qubit_slices...)
    A = similar(ψr)

    A[:, 1, :, 1, :] .= ψr[:, 1, :, 1, :]
    A[:, 2, :, 2, :] .= ψr[:, 2, :, 2, :]
    A[:, 1, :, 2, :] .= ψr[:, 2, :, 1, :]
    A[:, 2, :, 1, :] .= ψr[:, 1, :, 2, :]

    return reshape(A, :)
end


"""Tailored apply for a general ControlledGate"""
function _apply(cg::CircuitGate{M,ControlledGate{G}}, ψ::Vector{<:Complex}, N::Int) where {M,G}
    A = copy(ψ) # TODO: copy only what's needed
    A = reshape(A, fill(2, N)...)
    T = target_wires(cg.gate)
    itarget  = cg.iwire[1:T]
    icontrol = cg.iwire[T+1:M]

    new_iwire = Tuple(i - sum(icontrol .< i) for i in itarget)

    slice_target = Tuple(i in icontrol ? 2 : Colon() for i in 1:N)
    new_cg = CircuitGate(new_iwire, cg.gate.U)

    A[slice_target...] = reshape(apply(new_cg, reshape(A[slice_target...], :)),
                                 fill(2, N-M+T)...)
    return reshape(A, :)
end

"""
    apply(m::Moment, ψ::Vector{<:Complex}) 

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""
function apply(m::Moment, ψ::Vector{<:Complex}) 
    N_moment = size(m)
    l = length(ψ)::Int    
    N_ψ = Qaintessent.intlog2(l)
    l == 2^N_ψ || error("Vector length must be a power of 2")
    N_moment <= N_ψ || error("Moment affecting $N_moment qubits applied to $N_ψ qubits")
    _apply(m, ψ, N_ψ)   
end

function apply(m::Vector{Moment}, ψ::Vector{<:Complex})
    N_moment = maximum(size.(m))
    l = length(ψ)::Int    
    N_ψ = Qaintessent.intlog2(l)
    l == 2^N_ψ || error("Vector length must be a power of 2")
    N_moment <= N_ψ || error("Moment affecting $N_moment qubits applied to $N_ψ qubits")
    for moment in m
        ψ = _apply(moment, ψ, N_ψ)
    end
    ψ
end

function _apply(m::Moment, ψ::Vector{<:Complex}, N::Int)
    for gate in m
        ψ = _apply(gate, ψ, N)
    end
    return ψ
end

"""
    apply(m::MeasurementOperator, ψ::Vector{<:Complex}) 

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""

function apply(m::MeasurementOperator{M,G}, ψ::Vector{<:Complex}) where {M,G<:AbstractGate}
    N_moment = size(m)
    l = length(ψ)::Int    
    N_ψ = Qaintessent.intlog2(l)
    l == 2^N_ψ || error("Vector length must be a power of 2")
    N_moment <= N_ψ || error("MeasurementOperator affecting $N_moment qubits applied to $N_ψ qubits")
    c = circuit_gate((m.iwire...), m.operator)
    _apply(c, ψ, N_ψ)   
end

function apply(m::MeasurementOperator{M,G}, ψ::Vector{<:Complex}) where {M,G<:AbstractMatrix}
    N_moment = size(m)
    l = length(ψ)::Int    
    N_ψ = Qaintessent.intlog2(l)
    l == 2^N_ψ || error("Vector length must be a power of 2")
    N_moment <= N_ψ || error("MeasurementOperator affecting $N_moment qubits applied to $N_ψ qubits")
    apply(m.operator, ψ)   
end

# """
#     apply(cgc::CircuitGateChain{N}, ψ::Vector{<:Complex}) 

# Apply CircuitGateChain to quantum state vector.
# """
# function apply(cgc::CircuitGateChain{N}, ψ::AbstractVector{Complex}) 
#     creg = collect(Iterators.flatten(reverse.(cgc.creg)))
#     for moment in cgc.moments
#         for gate in moment
#             skip = false
#             for cntrl in gate.ccntrl
#                 if cntrl isa Expr
#                     if !eval(cntrl)
#                         skip = true
#                         break
#                     end
#                 else
#                     if creg[cntrl] != true
#                         skip = true
#                         break
#                     end
#                 end
#             end
#             if !skip
#                 ψ = apply(gate, ψ)
#             end
#         end
#     end
#     return ψ
# end


"""
    apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}

returns list of expectation values from measurement operators in `c.meas` after applying circuit gates in `c.cgc` on state vector of `N` qubits `ψ`
"""
function apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}
    length(ψ)::Int == 2^N || error("Size of vector `ψ` must match Circuit size of $(2^N)")
    ψl = deepcopy(ψ)
    for moment in c.moments
        ψl = _apply(moment, ψl, N)
    end
    ψr = apply.(c.meas, (ψl,)) 
    return real.(dot.((ψl,), ψr))
end

(c::Circuit)(ψ) = apply(c, ψ)


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
