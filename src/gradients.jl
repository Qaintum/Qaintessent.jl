
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

# TODO: derivatives of general RotationGate


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
