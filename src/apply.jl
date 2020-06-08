function bitswap(n::Int, p1::Int, p2::Int, N::Int)
    if p1 == p2
        return n
    end
    p1̄ = N-p1
    p2̄ = N-p2
    bit1 = (n >> p1̄) & 1
    bit2 = (n >> p2̄) & 1
    x = (bit1 ⊻ bit2)
    x = (x << p1̄) | (x << p2̄)
    return n ⊻ x
end

ror(x::Int, k::Int) = (x >>> (0x3f & k)) | (x << (0x3f & -k))
rol(x::Int, k::Int) = (x << (0x3f & k)) | (x >>> (0x3f & -k))

function swap(w::AbstractArray, p1::Int, p2::Int)
    tmp = w[p1]
    w[p1] = w[p2]
    w[p2] = tmp
    return w
end

"""
    apply(cg, ψ)

Apply general CircuitGate to quantum state vector
"""
function apply(cg::CircuitGate{M,N,G}, ψ::AbstractVector) where {M,N,G}
    W = length(cg.iwire)
    wires = [i for i in cg.iwire]

    gmat = matrix(cg.gate)
    # Use indices starting from 0 for easier calculation
    indices = 0:2^(N)-1

    for i in 1:W
        target = N-W+i
        indices = bitswap.(indices, wires[i], target, N)
        if target in wires
            target_index = Base.findfirst(x -> x==target, wires)
            wires = swap(wires, i, target_index)
        end
    end

    indices = sortperm(indices)

    ψ = ψ[indices]

    ψr = reshape(ψ, size(gmat)[1], :)

    ψr = gmat * ψr

    ψr = reshape(ψr, :)

    indices = sortperm(indices)
    return ψr[indices]
end


"""Tailored apply for XGate"""
function apply(cg::CircuitGate{1,N,XGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))

    return reshape(ψ[:,[2,1],:],:)
end

"""Tailored apply for YGate"""
function apply(cg::CircuitGate{1,N,YGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= -im.*ψ[:,2,:]
    A[:,2,:] .=  im.*ψ[:,1,:]
    return reshape(A,:)
end

"""Tailored apply for ZGate"""
function apply(cg::CircuitGate{1,N,ZGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= .- ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for HadamardGate"""
function apply(cg::CircuitGate{1,N,HadamardGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= (ψ[:,1,:] .+ ψ[:,2,:])./sqrt(2)
    A[:,2,:] .= (ψ[:,1,:] .- ψ[:,2,:])./sqrt(2)
    return reshape(A,:)
end

"""Tailored apply for SGate"""
function apply(cg::CircuitGate{1,N,SGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= im.*ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for SdagGate"""
function apply(cg::CircuitGate{1,N,SdagGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= -im.*ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for TGate"""
function apply(cg::CircuitGate{1,N,TGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(im*π/4).*ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for TdagGate"""
function apply(cg::CircuitGate{1,N,TdagGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(-im*π/4).*ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for RzGate"""
function apply(cg::CircuitGate{1,N,RzGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    θ = cg.gate.θ[]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= exp(-im*θ/2).*ψ[:,1,:]
    A[:,2,:] .= exp( im*θ/2).*ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for PhaseShiftGate"""
function apply(cg::CircuitGate{1,N,PhaseShiftGate}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(im*cg.gate.ϕ[]).*ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for a general single qubit gate"""
function apply(cg::CircuitGate{1,N,AbstractGate{1}}, ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    U = matrix(cg.gate)
    A[:,1,:] .= U[1,1] .* ψ[:,1,:] .+ U[1,2] .* ψ[:,2,:]
    A[:,2,:] .= U[2,1] .* ψ[:,1,:] .+ U[2,2] .* ψ[:,2,:]
    return reshape(A, :)
end

"""Tailored apply for SwapGate"""
function apply(cg::CircuitGate{2,N,SwapGate}, ψ::AbstractVector) where {N}

    i,j = cg.iwire
    i,j = i<j ? (i,j) : (j,i) #sort them
    qubit_slices = [i-1, 1, j-i-1, 1, N-j]

    ψr = reshape(ψ, 2 .^reverse(qubit_slices)...)
    A = similar(ψr)

    A[:,1,:,1,:] .= ψr[:,1,:,1,:]
    A[:,2,:,2,:] .= ψr[:,2,:,2,:]
    A[:,1,:,2,:] .= ψr[:,2,:,1,:]
    A[:,2,:,1,:] .= ψr[:,1,:,2,:]

    return reshape(A, :)
end

"""Tailored apply for a general ControlledGate"""
function apply(cg::CircuitGate{M,N,ControlledGate{T,M}}, ψ::AbstractVector) where {M,N,T}

    # M is the len of iwire, T is the number of target wires

    A = copy(ψ) # TODO: copy only what's needed
    A = reshape(A, fill(2, N)...)

    # Apply gate to qubits with control qubit = 2
    icontrol = cg.iwire[1:M-T]
    itarget =  cg.iwire[M-T+1:M]
    new_iwire= Tuple(i - sum(icontrol .< i) for i in itarget)
    slice_target = Tuple(i in icontrol ? 2 : Colon() for i in N:-1:1) # Why it has to go backwards?
    new_cg = CircuitGate(new_iwire, cg.gate.U, N-M+T)
    A[slice_target...] = reshape(apply(new_cg, reshape(A[slice_target...], :)),
                                 fill(2, N-M+T)...)
    return reshape(A, :)
end

"""
    apply(cgc, ψ)

Apply Moment to quantum state vector
"""
function apply(m::Moment{N}, ψ::AbstractVector) where {N}
    for gate in m
        ψ = apply(gate, ψ)
    end
    return ψ
end


"""
    apply(cgc, ψ)

Apply CircuitGateChain to quantum state vector.
"""
function apply(cgc::CircuitGateChain{N}, ψ::AbstractVector) where {N}
    for m in cgc.moments
        ψ = apply(m, ψ)
    end
    return ψ
end


"""
    apply(c, ψ)

Apply Circuit to quantum state vector
"""
function apply(c::Circuit{N}, ψ::AbstractVector) where {N}
    ψs = apply(c.cgc, ψ)
    return [real(dot(ψs, m*ψs)) for m in c.meas.mops]
end

(c::Circuit)(ψ) = apply(c, ψ)


"""Tailored apply to density matrix for XGate"""
function apply!(ρ::DensityMatrix{N}, cg::CircuitGate{1,N,XGate}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        # X I X =  I,
        # X X X =  X,
        # X Y X = -Y,
        # X Z X = -Z
        if pt[j] == 2 || pt[j] == 3
            ρ.v[i] *= -1
        end
    end
    return ρ
end

"""Tailored apply to density matrix for XGate"""
function apply(cg::CircuitGate{1,N,XGate}, ρ::DensityMatrix{N}) where {N}
    ρ = DensityMatrix{N}(copy(ρ.v))
    apply!(ρ, cg)
    return ρ
end


"""Tailored apply to density matrix for YGate"""
function apply!(ρ::DensityMatrix{N}, cg::CircuitGate{1,N,YGate}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        # Y I Y =  I,
        # Y X Y = -X,
        # Y Y Y =  Y,
        # Y Z Y = -Z
        if pt[j] == 1 || pt[j] == 3
            ρ.v[i] *= -1
        end
    end
    return ρ
end

"""Tailored apply to density matrix for YGate"""
function apply(cg::CircuitGate{1,N,YGate}, ρ::DensityMatrix{N}) where {N}
    ρ = DensityMatrix{N}(copy(ρ.v))
    apply!(ρ, cg)
    return ρ
end


"""Tailored apply to density matrix for ZGate"""
function apply!(ρ::DensityMatrix{N}, cg::CircuitGate{1,N,ZGate}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        # Z I Z =  I,
        # Z X Z = -X,
        # Z Y Z = -Y,
        # Z Z Z =  Z
        if pt[j] == 1 || pt[j] == 2
            ρ.v[i] *= -1
        end
    end
    return ρ
end

"""Tailored apply to density matrix for ZGate"""
function apply(cg::CircuitGate{1,N,ZGate}, ρ::DensityMatrix{N}) where {N}
    ρ = DensityMatrix{N}(copy(ρ.v))
    apply!(ρ, cg)
    return ρ
end


"""Tailored apply to density matrix for HadamardGate"""
function apply(cg::CircuitGate{1,N,HadamardGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    jbasefac = 4^(j-1)
    v = similar(ρ.v)
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        if pt[j] == 0
            # H I H = I
            v[i] = ρ.v[i]
        elseif pt[j] == 1
            # H Z H = X
            v[i] = ρ.v[i+2*jbasefac]
        elseif pt[j] == 2
            # H Y H = -Y
            v[i] = -ρ.v[i]
        elseif pt[j] == 3
            # H X H = Z
            v[i] = ρ.v[i-2*jbasefac]
        end
    end
    return DensityMatrix{N}(v)
end


"""Tailored apply to density matrix for SGate"""
function apply(cg::CircuitGate{1,N,SGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    jbasefac = 4^(j-1)
    v = similar(ρ.v)
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        if pt[j] == 0 || pt[j] == 3
            # S I S^† = I,
            # S Z S^† = Z
            v[i] = ρ.v[i]
        elseif pt[j] == 1
            # S Y S^† = -X
            v[i] = -ρ.v[i+jbasefac]
        elseif pt[j] == 2
            # S X S^† = Y
            v[i] = ρ.v[i-jbasefac]
        end
    end
    return DensityMatrix{N}(v)
end


"""Tailored apply to density matrix for SdagGate"""
function apply(cg::CircuitGate{1,N,SdagGate}, ρ::DensityMatrix{N}) where {N}
    # qubit index the gate acts on
    j = cg.iwire[1]
    jbasefac = 4^(j-1)
    v = similar(ρ.v)
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        if pt[j] == 0 || pt[j] == 3
            # S^† I S = I,
            # S^† Z S = Z
            v[i] = ρ.v[i]
        elseif pt[j] == 1
            # S^† Y S = X
            v[i] = ρ.v[i+jbasefac]
        elseif pt[j] == 2
            # S^† X S = -Y
            v[i] = -ρ.v[i-jbasefac]
        end
    end
    return DensityMatrix{N}(v)
end
