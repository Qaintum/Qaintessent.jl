function Qaintessent.apply(cg::CircuitGate{1,N,XGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))

    return reshape(ψ[:,[2,1],:],:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,YGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= -im.*ψ[:,2,:]
    A[:,2,:] .=  im.*ψ[:,1,:]
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,ZGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= -1 .*ψ[:,2,:]
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,HadamardGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= (ψ[:,1,:] .+ ψ[:,2,:])./sqrt(2)
    A[:,2,:] .= (ψ[:,1,:] .- ψ[:,2,:])./sqrt(2)
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,SGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= im.*ψ[:,2,:]
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,SdagGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= -im.*ψ[:,2,:]
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,TGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(im*π/4).*ψ[:,2,:]
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,TdagGate},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= ψ[:,1,:]
    A[:,2,:] .= exp(-im*π/4).*ψ[:,2,:]
    return reshape(A,:)
end

function Qaintessent.apply(cg::CircuitGate{1,N,AbstractGate{1}},ψ::AbstractVector) where {N}
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(N-i), 2, 2^(i-1))
    A = similar(ψ)
    A[:,1,:] .= U[1,1] .* ψ[:, 1, :] .+ U[1,2] .* ψ[:, 2, :]
    A[:,2,:] .= U[2,1] .* ψ[:, 1, :] .+ U[2,2] .* ψ[:, 2, :]
    return reshape(A,:)
end
