using Memoize

@memoize Dict function gradVector(N::Int)
    Statevector(N)
end

@memoize Dict function gradCircuit(c::Circuit)
    deepcopy(c)
end

# `backward` functions return gates storing derivatives

function backward!(dg::XGate, g::XGate, ::AbstractMatrix) end
function backward!(dg::YGate, g::YGate, ::AbstractMatrix) end
function backward!(dg::ZGate, g::ZGate, ::AbstractMatrix) end

function backward!(dg::HadamardGate, g::HadamardGate, ::AbstractMatrix) end


function backward!(dg::SGate, g::SGate, ::AbstractMatrix) end
function backward!(dg::TGate, g::TGate, ::AbstractMatrix) end

function backward!(dg::SdagGate, g::SdagGate, ::AbstractMatrix) end
function backward!(dg::TdagGate, g::TdagGate, ::AbstractMatrix) end

function backward!(dg::MatrixGate, g::MatrixGate, ::AbstractMatrix) end


function backward!(dg::RxGate, g::RxGate, Δ::AbstractMatrix)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    dg.θ .= sum(real([-s im * c; im * c -s] .* Δ))
    return
end

function backward!(dg::RyGate, g::RyGate, Δ::AbstractMatrix)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    dg.θ .= sum(real([-s -c; c -s] .* Δ))
    return
end

function backward!(dg::RzGate, g::RzGate, Δ::AbstractMatrix)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    eθ = exp(im * g.θ[] / 2)
    dg.θ .= real(im * eθ * Δ[1, 1] - im * conj(eθ) * Δ[2, 2])
    return
end


function backward!(dg::RotationGate, g::RotationGate, Δ::AbstractMatrix)
    θ = norm(g.nθ)
    if θ == 0
        for (i,p) in Enumerate([X,Y,Z])
            dg.nθ[i] = real(sum(conj(-im * p) .* Δ))
        end
    else
        n = g.nθ / θ
        c = cos(θ / 2)
        s = sin(θ / 2)
        dRθ = -s * I - im * c * pauli_vector(n...)
        dn = (I - reshape(kron(n, n), 3, 3)) / θ
        # using conjugated derivative matrix
        for i in 1:3
            dg.nθ[i] = 2 * real(sum(conj(0.5 * n[i] * dRθ - im * s * pauli_vector(dn[:, i]...)) .* Δ))
        end
    end
    return
end


function backward!(dg::PhaseShiftGate, g::PhaseShiftGate, Δ::AbstractMatrix)
    # using conjugated derivative matrix
    dg.ϕ .= 2 * real(-im * Base.exp(-im * g.ϕ[1]) * Δ[2, 2])
    return
end


function backward!(dg::SwapGate, g::SwapGate, ::AbstractMatrix) end


function backward!(dg::EntanglementXXGate, g::EntanglementXXGate, Δ::AbstractMatrix)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    dg.θ .= sum(real([-s 0 0 im * c; 0 -s im * c 0; 0 im * c -s 0; im * c 0 0 -s] .* Δ))
    return
end

function backward!(dg::EntanglementYYGate, g::EntanglementYYGate, Δ::AbstractMatrix)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    dg.θ .= sum(real([-s 0 0 -im * c; 0 -s im * c 0; 0 im * c -s 0; -im * c 0 0 -s] .* Δ))
    return
end

function backward!(dg::EntanglementZZGate, g::EntanglementZZGate, Δ::AbstractMatrix)
    eθ = exp(im * g.θ[] / 2)
    # using conjugated derivative matrix; factor 2 cancels 1/2 from θ/2
    dg.θ .= sum(real(im * eθ * Δ[1, 1] - im * conj(eθ) * Δ[2, 2] - im * conj(eθ) * Δ[3, 3] + im * eθ * Δ[4, 4]))
    return
end


function backward!(dg::ControlledGate{G}, g::ControlledGate{G}, Δ::AbstractMatrix) where {G}
    # Note: target qubits correspond to fastest varying indices
    backward!(dg.U, g.U, Δ[end-2^target_wires(g)+1:end, end-2^target_wires(g)+1:end])
    return
end

# `backward` functions return gates storing derivatives
function backward!(dcg::CircuitGate{M,G}, cg::CircuitGate{M,G}, ψ::Statevector, Δ::Statevector, N::Int) where {M,G}
    ρ = rdm(N, cg.iwire, Δ, ψ)
    backward!(dcg.gate, cg.gate, ρ)
    return
end

function backward!(dm::Moment, m::Moment, ψ::Statevector, Δ::Statevector, N::Int)
    i = length(m.gates)
    for cg in reverse(m.gates)
        Udag = Base.adjoint(cg)
        apply!(ψ, Udag)
        # backward step of quantum state
        backward!(dm[i], cg, ψ, Δ, N)
        apply!(Δ, Udag)
        i = i - 1
    end
    return Δ
end

function backward!(c::Circuit{N}, ψ::Statevector, Δ::Statevector) where {N}
    cl = gradCircuit(c)
    i = length(c.moments)
    for moment in reverse(c.moments)
        # backward step of quantum state
        Δ = backward!(cl[i], moment, ψ, Δ, N)
        i = i-1
    end
    return cl, Δ
end


"""
    gradients(c::Circuit{N}, ψ::AbstractVector, Δ::AbstractVector{<:Real}) where {N}

Perform a backward pass to compute gradients of a (fictitious) cost function with
respect to the circuit parameters of `c` and input wavefunction `ψ`. `Δ` contains
the cost function gradients with respect to the circuit outputs (measurement operator averages).
The gradients with respect to the circuit parameters are returned in a duplicate circuit;
the overall return value is the tuple (dc::Circuit{N}, dψ::AbstractVector).
"""
function gradients(c::Circuit{N}, ψ::Statevector, Δ::AbstractVector{<:Real}) where {N}
    # length of circuit output vector must match gradient vector
    @assert length(Δ) == length(c.meas)
    # forward pass through unitary gates
    ψl = gradVector(N)
    ψl .= ψ

    apply!(ψl, c.moments)
    # gradient (conjugated Wirtinger derivatives) of cost function with respect to ψ
    ψbar = Statevector(sum([Δ[i] * (sparse_matrix(c.meas[i]) * ψl) for i in 1:length(Δ)]))

    # # TODO: efficiently represent Kronecker product without explicitly storing matrix entries
    v = ψl.state
    dmeas = MeasurementOperator.([Δ[i] * reshape(kron(conj(v), v), length(v), length(v)) for i in 1:length(Δ)], (Tuple(1:N),))

    # # backward pass through unitary gates
    cl, ψbar = backward!(c, ψl, ψbar)
    cl.meas .= dmeas
    return cl, ψbar
end
