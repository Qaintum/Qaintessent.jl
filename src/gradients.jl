# derivatives of rotation gates with respect to rotation angle

struct dRxGate <: AbstractGate{1}
    θ::Real
end

function matrix(g::dRxGate)
    c = cos(g.θ/2)
    s = sin(g.θ/2)
    [-0.5*s -0.5im*c; -0.5im*c -0.5*s]
end

struct dRyGate <: AbstractGate{1}
    θ::Real
end

function matrix(g::dRyGate)
    c = cos(g.θ/2)
    s = sin(g.θ/2)
    [-0.5*s -0.5*c; 0.5*c -0.5*s]
end

struct dRzGate <: AbstractGate{1}
    θ::Real
end

function matrix(g::dRzGate)
    [-0.5im*exp(-im*g.θ/2) 0; 0 0.5im*exp(im*g.θ/2)]
end

dx(i, x, y, θ) = i == 1 ? 1/θ - x*y / θ^3 : - x*y / θ^3
dy(i, x, y, θ) = i == 2 ? - 1/θ + x*y / θ^3 : x*y / θ^3
dz(i, x, y, θ) = i == 3 ? 1/θ - x*y / θ^3 : - x*y / θ^3

function dpauli(i::Int, θ::Real, nθ::AbstractVector{<:Real})
    x = dx(i, nθ[1], nθ[i], θ)
    y = dy(i, nθ[2], nθ[i], θ)
    z = dz(i, nθ[3], nθ[i], θ)
    pauli_vector(x, y, z)
end

struct dRotationGateθ <: AbstractGate{1}
    θ::Real
    n::AbstractVector{<:Real}
end

struct dRotationGateG1 <: AbstractGate{1}
    θ::Real
    n::AbstractVector{<:Real}
end

struct dRotationGateG2 <: AbstractGate{1}
    θ::Real
    n::AbstractVector{<:Real}
end

struct dRotationGateG3 <: AbstractGate{1}
    θ::Real
    n::AbstractVector{<:Real}
end

struct dRotationGateX <: AbstractGate{1}
    nθ::AbstractVector{<:Real}
end

struct dRotationGateY <: AbstractGate{1}
    nθ::AbstractVector{<:Real}
end

struct dRotationGateZ <: AbstractGate{1}
    nθ::AbstractVector{<:Real}
end

function matrix(g::dRotationGateθ)

    -0.5*sin(g.θ/2)*I - 0.5im*cos(g.θ/2) * conj(pauli_vector(g.n...))
end

function matrix(g::dRotationGateG1)

    -im*sin(g.θ/2) * pauli_vector(1, 0, 0)
end

function matrix(g::dRotationGateG2)

    -im*sin(g.θ/2) * pauli_vector(0, -1, 0)
end

function matrix(g::dRotationGateG3)

    -im*sin(g.θ/2) * pauli_vector(0, 0, 1)
end

function matrix(g::dRotationGateX)
    θ = norm(g.nθ)
    n = g.nθ/θ

    - g.nθ[1]/2θ*sin(θ/2)*I - im*g.nθ[1]/2θ*cos(θ/2) * conj(pauli_vector(n...)) - im*sin(θ/2) * dpauli(1, θ, g.nθ)
end

function matrix(g::dRotationGateY)
    θ = norm(g.nθ)
    n = g.nθ/θ

    - g.nθ[2]/2θ*sin(θ/2)*I - im*g.nθ[2]/2θ*cos(θ/2) * conj(pauli_vector(n...)) - im*sin(θ/2) * dpauli(2, θ, g.nθ)
end

function matrix(g::dRotationGateZ)
    θ = norm(g.nθ)
    n = g.nθ/θ

    - g.nθ[3]/2θ*sin(θ/2)*I - im*g.nθ[3]/2θ*cos(θ/2) * conj(pauli_vector(n...)) - im*sin(θ/2) * dpauli(3, θ, g.nθ)
end


struct dPhaseShiftGate <: AbstractGate{1}
    ϕ::Real
end

matrix(g::dPhaseShiftGate) = [0 0; 0 im*exp(im*g.ϕ)]


struct dControlledGate{M,N} <: AbstractGate{N}
    dU::AbstractGate{M}
end

function matrix(g::dControlledGate{M,N}) where {M,N}
    dUmat = matrix(g.dU)
    dCU = sparse([], [], eltype(dUmat)[], 2^N, 2^N)
    # Note: following the ordering convention of `kron` here, i.e.,
    # target qubits correspond to fastest varying index
    dCU[end-size(dUmat,1)+1:end, end-size(dUmat,2)+1:end] = dUmat
    return dCU
end


# collect first derivatives

function derivatives(g::RxGate)
    [dRxGate(g.θ)]
end

function derivatives(g::RyGate)
    [dRyGate(g.θ)]
end

function derivatives(g::RzGate)
    [dRzGate(g.θ)]
end

function derivatives(g::PhaseShiftGate)
    [dPhaseShiftGate(g.ϕ)]
end

function derivatives(g::RotationGate)
    # [dRotationGateθ(g.θ, g.n), dRotationGateG1(g.θ, g.n), dRotationGateG2(g.θ, g.n), dRotationGateG3(g.θ, g.n)]
    [dRotationGateX(g.nθ), dRotationGateY(g.nθ), dRotationGateZ(g.nθ)]
end

function derivatives(g::ControlledGate{M,N}) where {M,N}
    [dControlledGate{M,N}(dU) for dU in derivatives(g.U)]
end


# default: no derivatives
derivatives(g::AbstractGate{N}) where {N} = []

function derivatives(cg::CircuitGate{M,N}) where {M,N}
    [CircuitGate{M,N}(cg.iwire, dg) for dg in derivatives(cg.gate)]
end


function gradients(c::Circuit{N}, ψ::AbstractVector, Δ::AbstractVector{<:Real}) where {N}
    # length of circuit output vector must match gradient vector
    @assert length(Δ) == length(c.meas.mops)

    # forward pass through unitary gates
    ψ = apply(c.cgc, ψ)
    # gradient (conjugated Wirtinger derivatives) of cost function with respect to ψ
    ψbar = sum([Δ[i] * (c.meas.mops[i]*ψ) for i in 1:length(Δ)])

    # assuming for now that all parameters are real-valued
    grads = Real[]
    for cg in reverse(c.cgc.gates)
        Udag = Base.adjoint(cg)
        # backward step of quantum state
        ψ = apply(Udag, ψ)
        # use reverse here to counteract final reverse
        for dcg in reverse(derivatives(cg))
            push!(grads, 2*real(dot(ψ, conj(matrix(dcg))*ψbar)))
        end
        # backward step of quantum state gradient
        ψbar = apply(Udag, ψbar)
    end
    return reverse(grads)
end
