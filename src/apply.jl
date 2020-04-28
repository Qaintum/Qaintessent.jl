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

Apply CircuitGateChain to quantum state vector
"""
function apply(cgc::CircuitGateChain{N}, ψ::AbstractVector) where {N}
    for moment in cgc.moments
        for gate in moment
            ψ = apply(gate, ψ)
        end
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
