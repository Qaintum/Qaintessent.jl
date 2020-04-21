"""Tailored apply for XGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,XGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))

    return reshape(ψ[:,[2,1],:],:)
end

"""Tailored apply for YGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,YGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= -im.*ψ[:,2,:]
    A[:,2,:] .=  im.*ψ[:,1,:]
    return reshape(A,:)
end

"""Tailored apply for ZGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,ZGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= -1 .*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for HadamardGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,HadamardGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= (ψ[:,1,:] .+ ψ[:,2,:])./sqrt(2)
    A[:,2,:] .= (ψ[:,1,:] .- ψ[:,2,:])./sqrt(2)
    return reshape(A,:)
end

"""Tailored apply for SGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,SGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= im.*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for SdagGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,SdagGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= -im.*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for TGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,TGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(im*π/4).*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for TdagGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,TdagGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(-im*π/4).*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for RzGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,RzGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    θ = cg.gate.θ[]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= exp(-im*θ/2).*ψ[:,1,:]
    A[:,2,:] .= exp( im*θ/2).*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for PhaseShiftGate"""
function Qaintessent.apply(cg::CircuitGate{1,N,PhaseShiftGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp( im*cg.gate.ϕ[]).*ψ[:,2,:]
    return reshape(A,:)
end

"""Tailored apply for a general single qubit gate"""
function Qaintessent.apply(cg::CircuitGate{1,N,AbstractGate{1}},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= U[1,1] .* ψ[:, 1, :] .+ U[1,2] .* ψ[:, 2, :]
    A[:,2,:] .= U[2,1] .* ψ[:, 1, :] .+ U[2,2] .* ψ[:, 2, :]
    return reshape(A,:)
end

"""Tailored apply for SwapGate"""
function Qaintessent.apply(cg::CircuitGate{2,N,SwapGate}, ψ::AbstractVector) where {N}

    i,j = cg.iwire
    i,j = i<j ? (i,j) : (j,i) #sort them
    qubit_slices = [i-1, 1, j-i-1, 1, N-j]

    ψr = reshape(ψ, 2 .^reverse(qubit_slices)...)
    A = similar(ψr)

    A[:,1,:,1,:] .= ψr[:,1,:,1,:]
    A[:,2,:,2,:] .= ψr[:,2,:,2,:]
    A[:,1,:,2,:] .= ψr[:,2,:,1,:]
    A[:,2,:,1,:] .= ψr[:,1,:,2,:]

    return reshape(A,:)
end

"""Tailored apply for a general ControlledGate"""
function Qaintessent.apply(cg::CircuitGate{M,N,ControlledGate{K,M}}, ψ::AbstractVector) where {M,N,K}

    # M is the len of iwire, K is the number of target wires

    A = copy(ψ) # TODO: copy only what's needed
    A = reshape(A, fill(2,N)...)

    # Apply gate to qubits with control qubit = 2
    icontrol = cg.iwire[1:M-K]
    itarget =  cg.iwire[M-K+1:M]
    new_iwire= Tuple(i - sum(icontrol .< i) for i in itarget)
    slice_target = Tuple(i in icontrol ? 2 : Colon() for i in N:-1:1) # Why it has to go backwards?
    new_cg = CircuitGate{K,N-M+K}(new_iwire, cg.gate.U)
    A[slice_target...] = reshape(apply(new_cg, reshape(A[slice_target...],:)),
                                 fill(2,N-M+K)...)

    return reshape(A,:)
end
