
# `backward` functions return gates storing derivatives

backward(g::XGate, ::Matrix{ComplexQ}) = g
backward(g::YGate, ::Matrix{ComplexQ}) = g
backward(g::ZGate, ::Matrix{ComplexQ}) = g


backward(g::HadamardGate, ::Matrix{ComplexQ}) = g


backward(g::SGate, ::Matrix{ComplexQ}) = g
backward(g::TGate, ::Matrix{ComplexQ}) = g

backward(g::SdagGate, ::Matrix{ComplexQ}) = g
backward(g::TdagGate, ::Matrix{ComplexQ}) = g

backward(g::MatrixGate, ::Matrix{ComplexQ}) = g


function backward(g::RxGate, Δ::Matrix{ComplexQ})
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    RxGate(sum(real([-s im * c; im * c -s] .* Δ)))
end

function backward(g::RyGate, Δ::Matrix{ComplexQ})
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    RyGate(sum(real(ComplexQ[-s -c; c -s] .* Δ)))
end

function backward(g::RzGate, Δ::Matrix{ComplexQ})
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    eθ = exp(im * g.θ[] / 2)
    RzGate(real(im * eθ * Δ[1, 1] - im * conj(eθ) * Δ[2, 2]))
end


function backward(g::RotationGate, Δ::Matrix{ComplexQ})
    θ = norm(g.nθ)
    if θ == 0
        σ = Matrix{ComplexF64}[matrix(X), matrix(Y), matrix(Z)]
        RotationGate([real(sum(conj(-im * σ[i]) .* Δ)) for i in 1:3])
    else
        n = g.nθ / θ
        c = cos(θ / 2)
        s = sin(θ / 2)
        dRθ = -s * I - im * c * pauli_vector(n...)
        dn = (I - reshape(kron(n, n), 3, 3)) / θ
        # using conjugated derivative matrix
        RotationGate([2 * real(sum(conj(0.5 * n[i] * dRθ - im * s * pauli_vector(dn[:, i]...)) .* Δ)) for i in 1:3])
    end
end


function backward(g::PhaseShiftGate, Δ::Matrix{ComplexQ})
    # using conjugated derivative matrix
    PhaseShiftGate(2 * real(-im * Base.exp(-im * g.ϕ[1]) * Δ[2, 2]))
end


backward(g::SwapGate, ::Matrix{ComplexQ}) = g


function backward(g::EntanglementXXGate, Δ::Matrix{ComplexQ})
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    EntanglementXXGate(sum(real([-s 0 0 im * c; 0 -s im * c 0; 0 im * c -s 0; im * c 0 0 -s] .* Δ)))
end

function backward(g::EntanglementYYGate, Δ::Matrix{ComplexQ})
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    EntanglementYYGate(sum(real([-s 0 0 -im * c; 0 -s im * c 0; 0 im * c -s 0; -im * c 0 0 -s] .* Δ)))
end

function backward(g::EntanglementZZGate, Δ::Matrix{ComplexQ})
    eθ = exp(im * g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    EntanglementZZGate(sum(real(im * eθ * Δ[1, 1] - im * conj(eθ) * Δ[2, 2] - im * conj(eθ) * Δ[3, 3] + im * eθ * Δ[4, 4])))
end


function backward(g::ControlledGate{G}, Δ::Matrix{ComplexQ}) where {G}
    # Note: target qubits correspond to fastest varying indices
    ControlledGate(backward(g.U, Δ[end-2^target_wires(g)+1:end, end-2^target_wires(g)+1:end]), control_wires(g))
end


function backward(cg::CircuitGate{M,G}, ψ::Vector{T}, Δ::Vector{T}, N::Int) where {M,G,T}
    ρ = rdm(N, cg.iwire, Δ, ψ)
    CircuitGate{M,G}(cg.iwire, backward(cg.gate, ρ))
end

function backward(m::Moment, ψ::Vector{T}, Δ::Vector{T}, N::Int) where {T}
    gates = AbstractCircuitGate[]
    for cg in reverse(m.gates)
        Udag = Base.adjoint(cg)
        ψ = apply(ψ, Udag)
        # backward step of quantum state
        pushfirst!(gates, backward(cg, ψ, Δ, N))
        Δ = apply(Δ, Udag)
    end
    return Moment(gates), ψ, Δ
end

function backward(moments::Vector{Moment}, ψ::Vector{T}, Δ::Vector{T}, N::Int) where {T}
    dmoments = Moment[]
    for moment in reverse(moments)
        # backward step of quantum state
        (m, ψ, Δ) = backward(moment, ψ, Δ, N)
        pushfirst!(dmoments, m)
    end
    return dmoments, Δ
end


"""
    gradients(c::Circuit{N}, ψ::AbstractVector, Δ::AbstractVector{<:Real}) where {N}

Perform a backward pass to compute gradients of a (fictitious) cost function with
respect to the circuit parameters of `c` and input wavefunction `ψ`. `Δ` contains
the cost function gradients with respect to the circuit outputs (measurement operator averages).
The gradients with respect to the circuit parameters are returned in a duplicate circuit;
the overall return value is the tuple (dc::Circuit{N}, dψ::AbstractVector).
"""
function gradients(c::Circuit{N}, ψ::Vector{ComplexQ}, Δ::AbstractArray) where {N}
    # length of circuit output vector must match gradient vector
    @assert length(Δ) == length(c.meas)
    # forward pass through unitary gates
    ψ = apply(ψ, c.moments)
    # gradient (conjugated Wirtinger derivatives) of cost function with respect to ψ
    ψbar = sum([Δ[i] * (sparse_matrix(c.meas[i]) * ψ) for i in 1:length(Δ)])
    # backward pass through unitary gates
    dcgc, ψbar = backward(c.moments, ψ, ψbar, N)
    # TODO: efficiently represent Kronecker product without explicitly storing matrix entries
    dmeas = MeasurementOperator.([Δ[i] * reshape(kron(conj(ψ), ψ), length(ψ), length(ψ)) for i in 1:length(Δ)], (Tuple(1:N),))
    return Circuit{N}(dcgc, dmeas), ψbar
end
