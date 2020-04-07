

# `backward` functions return an `IdDict` with the corresponding variable hash as key

function backward(g::RxGate, Δ::AbstractMatrix)
    c = cos(g.θ[1]/2)
    s = sin(g.θ[1]/2)
    # using conjugated derivative matrix
    IdDict(g.θ => [2*sum(real([-0.5*s 0.5im*c; 0.5im*c -0.5*s] .* Δ))])
end

function backward(g::RyGate, Δ::AbstractMatrix)
    c = cos(g.θ[1]/2)
    s = sin(g.θ[1]/2)
    IdDict(g.θ => [2*sum(real([-0.5*s -0.5*c; 0.5*c -0.5*s] .* Δ))])
end

function backward(g::RzGate, Δ::AbstractMatrix)
    # using conjugated derivative matrix
    IdDict(g.θ => [2*real(0.5im*exp(im*g.θ[1]/2)*Δ[1, 1] - 0.5im*exp(-im*g.θ[1]/2)*Δ[2, 2])])
end

dij(i,j, x, y, θ) = i == j  ? 1/θ - x*y / θ^3 : - x*y / θ^3

function dpauli(i::Int, θ::Real, nθ::AbstractVector{<:Real})
    v = [ dij(i, j, nθ[j], nθ[i], θ) for j in 1:3 ]
    pauli_vector(v...)
end

function backward(g::RotationGate, Δ::AbstractMatrix)
    θ = norm(g.nθ)
    n = g.nθ/θ

    p = pauli_vector(n...)
    dp = [dpauli(i, θ, g.nθ) for i in 1:3]

    mat1 = - g.nθ[1]/2θ*sin(θ/2)*I - im*g.nθ[1]/2θ*cos(θ/2) * p - im*sin(θ/2) * dp[1]
    mat2 = - g.nθ[2]/2θ*sin(θ/2)*I - im*g.nθ[2]/2θ*cos(θ/2) * p - im*sin(θ/2) * dp[2]
    mat3 = - g.nθ[3]/2θ*sin(θ/2)*I - im*g.nθ[3]/2θ*cos(θ/2) * p - im*sin(θ/2) * dp[3]
    # using conjugated derivative matrix
    IdDict(
    g.nθ => [2*real(sum( conj(mat1) .* Δ )),
             2*real(sum( conj(mat2) .* Δ )),
             2*real(sum( conj(mat3) .* Δ ))]
    )
end


function backward(g::PhaseShiftGate, Δ::AbstractMatrix)
    # using conjugated derivative matrix
    IdDict(g.ϕ => [2*real(-im*exp(-im*g.ϕ[1])*Δ[2, 2])])
end


function backward(g::ControlledGate{M,N}, Δ::AbstractMatrix) where {M,N}
    # Note: following the ordering convention of `kron` here, i.e.,
    # target qubits correspond to fastest varying index
    backward(g.U, Δ[end-2^M+1:end, end-2^M+1:end])
end


# default: no gradients
backward(g::AbstractGate{N}, Δ::AbstractMatrix) where {N} = IdDict()


function backward(cg::CircuitGate{M,N}, ψ::AbstractVector, Δ::AbstractVector) where {M,N}
    ρ = rdm(N, cg.iwire, Δ, ψ)
    backward(cg.gate, ρ)
end


function backward(cgc::CircuitGateChain{N}, ψ::AbstractVector, Δ::AbstractVector) where {N}
    grads = IdDict()
    for cg in reverse(cgc.gates)
        Udag = Base.adjoint(cg)
        # backward step of quantum state
        ψ = apply(Udag, ψ)
        grads = merge(backward(cg, ψ, Δ), grads)
        # backward step of quantum state gradient
        Δ = apply(Udag, Δ)
    end
    return grads, Δ
end


function gradients(c::Circuit{N}, ψ::AbstractVector, Δ::AbstractVector{<:Real}) where {N}
    # length of circuit output vector must match gradient vector
    @assert length(Δ) == length(c.meas.mops)
    # forward pass through unitary gates
    ψ = apply(c.cgc, ψ)
    # gradient (conjugated Wirtinger derivatives) of cost function with respect to ψ
    ψbar = sum([Δ[i] * (c.meas.mops[i]*ψ) for i in 1:length(Δ)])
    return backward(c.cgc, ψ, ψbar)[1]
end
