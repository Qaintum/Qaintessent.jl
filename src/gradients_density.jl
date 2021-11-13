# `backward_density` functions return gates storing derivatives

# The logic for gradient computation here is analogous to `gradients.jl`
# when interpreting the coefficient vector of DensityMatrix as statevector.


backward_density(g::XGate, ::AbstractMatrix) = g
backward_density(g::YGate, ::AbstractMatrix) = g
backward_density(g::ZGate, ::AbstractMatrix) = g


backward_density(g::HadamardGate, ::AbstractMatrix) = g


backward_density(g::SGate, ::AbstractMatrix) = g
backward_density(g::TGate, ::AbstractMatrix) = g

backward_density(g::SdagGate, ::AbstractMatrix) = g
backward_density(g::TdagGate, ::AbstractMatrix) = g


function backward_density(g::RxGate, Δ::AbstractMatrix)
    c = cos(g.θ[])
    s = sin(g.θ[])
    RxGate(-s*(Δ[3, 3] + Δ[4, 4]) + c*(Δ[4, 3] - Δ[3, 4]))
end

function backward_density_mixed_add(g::RxGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    RxGate(0.5*(-s*tr(Δ) + c*(Δ[4, 3] - Δ[3, 4])))
end

function backward_density_mixed_sub(g::RxGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    RxGate(0.5*c*(Δ[1, 2] + Δ[2, 1]))
end


function backward_density(g::RyGate, Δ::AbstractMatrix)
    c = cos(g.θ[])
    s = sin(g.θ[])
    RyGate(-s*(Δ[2, 2] + Δ[4, 4]) + c*(Δ[2, 4] - Δ[4, 2]))
end

function backward_density_mixed_add(g::RyGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    RyGate(0.5*(-s*tr(Δ) + c*(Δ[2, 4] - Δ[4, 2])))
end

function backward_density_mixed_sub(g::RyGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    RyGate(0.5*c*(Δ[1, 3] + Δ[3, 1]))
end


function backward_density(g::RzGate, Δ::AbstractMatrix)
    c = cos(g.θ[])
    s = sin(g.θ[])
    RzGate(-s*(Δ[2, 2] + Δ[3, 3]) + c*(Δ[3, 2] - Δ[2, 3]))
end

function backward_density_mixed_add(g::RzGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    RzGate(0.5*(-s*tr(Δ) + c*(Δ[3, 2] - Δ[2, 3])))
end

function backward_density_mixed_sub(g::RzGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    RzGate(0.5*c*(Δ[1, 4] + Δ[4, 1]))
end


function backward_density(g::RotationGate, Δ::AbstractMatrix)
    θ = norm(g.nθ)
    if θ == 0
        RotationGate([Δ[4, 3] - Δ[3, 4], Δ[2, 4] - Δ[4, 2], Δ[3, 2] - Δ[2, 3]])
    else
        n = g.nθ/θ
        c = cos(θ)
        s = sin(θ)
        # cross-product matrix
        K = [0 -n[3] n[2]; n[3] 0 -n[1]; -n[2] n[1] 0]
        # -K^2
        K2 = I - reshape(kron(n, n), 3, 3)
        dRθ = c*K - s*K2
        dn = K2 / θ
        dRn1 = [2*(1-c)*n[1] (1-c)*n[2] (1-c)*n[3]; (1-c)*n[2] 0 -s; (1-c)*n[3] s 0]
        dRn2 = [0 (1-c)*n[1] s; (1-c)*n[1] 2*(1-c)*n[2] (1-c)*n[3]; -s (1-c)*n[3] 0]
        dRn3 = [0 -s (1-c)*n[1]; s 0 (1-c)*n[2]; (1-c)*n[1] (1-c)*n[2] 2*(1-c)*n[3]]
        RotationGate([sum((n[i]*dRθ + (dRn1*dn[1, i] + dRn2*dn[2, i] + dRn3*dn[3, i])) .* Δ[2:4, 2:4]) for i in 1:3])
    end
end

function backward_density_mixed_add(g::RotationGate, Δ::AbstractMatrix)
    dsn = [Δ[4, 3] - Δ[3, 4], Δ[2, 4] - Δ[4, 2], Δ[3, 2] - Δ[2, 3]]
    θ = norm(g.nθ)
    if θ == 0
        return RotationGate(0.5 * dsn)
    else
        n = g.nθ/θ
        c = cos(0.5*θ)
        s = sin(0.5*θ)
        dnθ = -0.5*s * tr(Δ) * n + (0.5*c - s/θ) * dot(n, dsn) * n + (s/θ) * dsn
        return RotationGate(dnθ)
    end
end

function backward_density_mixed_sub(g::RotationGate, Δ::AbstractMatrix)
    dsn = [Δ[1, 2] + Δ[2, 1], Δ[1, 3] + Δ[3, 1], Δ[1, 4] + Δ[4, 1]]
    θ = norm(g.nθ)
    if θ == 0
        return RotationGate(0.5 * dsn)
    else
        n = g.nθ/θ
        c = cos(0.5*θ)
        s = sin(0.5*θ)
        dnθ = (0.5*c - s/θ) * dot(n, dsn) * n + (s/θ) * dsn
        return RotationGate(dnθ)
    end
end


function backward_density(g::PhaseShiftGate, Δ::AbstractMatrix)
    # agrees with backward pass for Rz(θ) since global prefactor cancels
    c = cos(g.ϕ[])
    s = sin(g.ϕ[])
    PhaseShiftGate(sum([-s -c; c -s] .* Δ[2:3, 2:3]))
end

function backward_density_mixed_add(g::PhaseShiftGate, Δ::AbstractMatrix)
    c = cos(g.ϕ[])
    s = sin(g.ϕ[])
    PhaseShiftGate(0.5*(s*((Δ[1, 4] + Δ[4, 1]) - tr(Δ)) + c*(Δ[3, 2] - Δ[2, 3])))
end

function backward_density_mixed_sub(g::PhaseShiftGate, Δ::AbstractMatrix)
    c = cos(g.ϕ[])
    s = sin(g.ϕ[])
    PhaseShiftGate(0.5*(c*((Δ[1, 4] + Δ[4, 1]) - tr(Δ)) + s*(Δ[2, 3] - Δ[3, 2])))
end


backward_density(g::SwapGate, Δ::AbstractMatrix) = g


function backward_density(g::EntanglementXXGate, Δ::AbstractMatrix)
    c = cos(g.θ[])
    s = sin(g.θ[])
    EntanglementXXGate(-s*(Δ[3, 3] + Δ[8, 8] + Δ[4, 4] + Δ[7, 7] + Δ[9,  9] + Δ[14, 14] + Δ[10, 10] + Δ[13, 13])
                       +c*(Δ[8, 3] - Δ[3, 8] + Δ[4, 7] - Δ[7, 4] + Δ[14, 9] - Δ[ 9, 14] + Δ[13, 10] - Δ[10, 13]))
end

function backward_density_mixed_add(g::EntanglementXXGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    EntanglementXXGate(0.5*(-s*tr(Δ) + c*(Δ[4, 7] - Δ[7, 4] + Δ[8, 3] - Δ[3, 8] + Δ[14, 9] - Δ[9, 14] + Δ[13, 10] - Δ[10, 13])))
end

function backward_density_mixed_sub(g::EntanglementXXGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    EntanglementXXGate(0.5*c*(Δ[6, 1] + Δ[1, 6] + Δ[5, 2] + Δ[2, 5] + Δ[15, 12] + Δ[12, 15] - Δ[16, 11] - Δ[11, 16]))
end


function backward_density(g::EntanglementYYGate, Δ::AbstractMatrix)
    c = cos(g.θ[])
    s = sin(g.θ[])
    EntanglementYYGate(-s*(Δ[2,  2] + Δ[12, 12] + Δ[10, 10] + Δ[4,  4] + Δ[5,  5] + Δ[15, 15] + Δ[7,  7] + Δ[13, 13])
                       +c*(Δ[2, 12] - Δ[12,  2] + Δ[10,  4] - Δ[4, 10] + Δ[5, 15] - Δ[15,  5] + Δ[7, 13] - Δ[13,  7]))
end

function backward_density_mixed_add(g::EntanglementYYGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    EntanglementYYGate(0.5*(-s*tr(Δ) + c*(Δ[10, 4] - Δ[4, 10] + Δ[2, 12] - Δ[12, 2] + Δ[5, 15] - Δ[15, 5] + Δ[7, 13] - Δ[13, 7])))
end

function backward_density_mixed_sub(g::EntanglementYYGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    EntanglementYYGate(0.5*c*(Δ[11, 1] + Δ[1, 11] + Δ[9, 3] + Δ[3, 9] + Δ[14, 8] + Δ[8, 14] - Δ[16, 6] - Δ[6, 16]))
end


function backward_density(g::EntanglementZZGate, Δ::AbstractMatrix)
    c = cos(g.θ[])
    s = sin(g.θ[])
    EntanglementZZGate(-s*(Δ[3,  3] + Δ[14, 14] + Δ[15, 15] + Δ[2,  2] + Δ[12, 12] + Δ[5,  5] + Δ[9, 9] + Δ[8, 8])
                       +c*(Δ[3, 14] - Δ[14,  3] + Δ[15,  2] - Δ[2, 15] + Δ[12,  5] - Δ[5, 12] + Δ[9, 8] - Δ[8, 9]))
end

function backward_density_mixed_add(g::EntanglementZZGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    EntanglementZZGate(0.5*(-s*tr(Δ) + c*(Δ[15, 2] - Δ[2, 15] + Δ[3, 14] - Δ[14, 3] + Δ[12, 5] - Δ[5, 12] + Δ[9, 8] - Δ[8, 9])))
end

function backward_density_mixed_sub(g::EntanglementZZGate, Δ::AbstractMatrix)
    c = cos(0.5*g.θ[])
    s = sin(0.5*g.θ[])
    EntanglementZZGate(0.5*c*(Δ[16, 1] + Δ[1, 16] + Δ[13, 4] + Δ[4, 13] + Δ[10, 7] + Δ[7, 10] - Δ[6, 11] - Δ[11, 6]))
end


# shortcuts for non-parametric gates

backward_density(g::ControlledGate{XGate}, ::AbstractMatrix) = g
backward_density(g::ControlledGate{YGate}, ::AbstractMatrix) = g
backward_density(g::ControlledGate{ZGate}, ::AbstractMatrix) = g

backward_density(g::ControlledGate{HadamardGate}, ::AbstractMatrix) = g

backward_density(g::ControlledGate{SGate}, ::AbstractMatrix) = g
backward_density(g::ControlledGate{TGate}, ::AbstractMatrix) = g
backward_density(g::ControlledGate{SdagGate}, ::AbstractMatrix) = g
backward_density(g::ControlledGate{TdagGate}, ::AbstractMatrix) = g

backward_density(g::ControlledGate{SwapGate}, ::AbstractMatrix) = g


function backward_density(g::ControlledGate{G}, Δ::AbstractMatrix) where {G}

    # number of target wires
    T = target_wires(g)

    # number of control wires
    C = control_wires(g)

    # for a single control qubit:
    # controlled-U = |0><0| x I + |1><1| x U = I + |1><1| x (U - I)
    # (generalizes to several control qubits)

    # mixed term |1><1| x (U - I) ρ + ρ |1><1| x (U† - I)
    Δadd = zeros(eltype(Δ), 4^T, 4^T)
    Δsub = zeros(eltype(Δ), 4^T, 4^T)
    eo = binary_digits(C, 0)
    for p in 0:2^C-1
        binary_digits!(eo, p)
        Δr = reshape(Δ, 4^T, 4^C, 4^T, 4^C)
        for j in 1:C
            # use another variable for type stability
            Δt = reshape(Δr, 4^T, 4, 4^(C-j), 4^T, 4, 4^(C-j))
            if eo[j] == 0
                Δr = 0.5 * (Base.view(Δt, :, 1, :, :, 1, :) + Base.view(Δt, :, 2, :, :, 2, :) + Base.view(Δt, :, 3, :, :, 3, :) + Base.view(Δt, :, 4, :, :, 4, :) - Base.view(Δt, :, 1, :, :, 4, :) - Base.view(Δt, :, 4, :, :, 1, :))
            else
                Δr = 0.5 * (Base.view(Δt, :, 3, :, :, 2, :) - Base.view(Δt, :, 2, :, :, 3, :))
            end
        end
        nodd::Int = sum(eo)
        # sign factor from even powers of i
        signfac = 1 - 2*(((nodd + 1) ÷ 2) % 2)
        if iseven(nodd)
            Δadd += (2*signfac) * reshape(Δr, 4^T, 4^T)
        else
            Δsub += (2*signfac) * reshape(Δr, 4^T, 4^T)
        end
    end
    # backward_density_mixed_add(g.U, Δadd) delayed until later
    dU = backward_density_mixed_sub(g.U, Δsub)

    # conjugation by |1><1| x (U - I)
    Δr = reshape(Δ, 4^T, 4^C, 4^T, 4^C)
    for j in 1:C
        # use another variable for type stability
        Δt = reshape(Δr, 4^T, 4, 4^(C-j), 4^T, 4, 4^(C-j))
        Δr = 0.5 * (Base.view(Δt, :, 1, :, :, 1, :) + Base.view(Δt, :, 4, :, :, 4, :) - Base.view(Δt, :, 1, :, :, 4, :) - Base.view(Δt, :, 4, :, :, 1, :))
    end
    dU += backward_density(g.U, reshape(Δr, 4^T, 4^T))
    Δadd -= 2 * reshape(Δr, 4^T, 4^T)
    dU += backward_density_mixed_add(g.U, Δadd)

    return ControlledGate{G}(dU, g.M)
end


backward_density(g::MatrixGate, ::AbstractMatrix) = g


# shortcuts for non-parametric gates

backward_density(cg::CircuitGate{M,XGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg
backward_density(cg::CircuitGate{M,YGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg
backward_density(cg::CircuitGate{M,ZGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg

backward_density(cg::CircuitGate{M,HadamardGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg

backward_density(cg::CircuitGate{M,SGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg
backward_density(cg::CircuitGate{M,TGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg
backward_density(cg::CircuitGate{M,SdagGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg
backward_density(cg::CircuitGate{M,TdagGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg

backward_density(cg::CircuitGate{M,SwapGate}, ::DensityMatrix, ::DensityMatrix) where {M} = cg

function backward_density(cg::CircuitGate{M,G}, ρ::DensityMatrix, Δ::DensityMatrix) where {M,G}
    @assert ρ.N == Δ.N
    CircuitGate{M,G}(cg.iwire, backward_density(cg.gate, rdm(ρ.N, cg.iwire, Δ.v, ρ.v, 4)))
end


function backward_density(m::Moment, ρ::DensityMatrix, Δ::DensityMatrix)
    dgates = CircuitGate[]
    for cg in reverse(m.gates)
        Udag = Base.adjoint(cg)
        apply!(ρ, Udag)
        # backward step of quantum state
        pushfirst!(dgates, backward_density(cg, ρ, Δ))
        apply!(Δ, Udag)
    end
    return Moment(dgates), ρ, Δ
end

function backward_density(moments::Vector{Moment}, ρ::DensityMatrix, Δ::DensityMatrix)
    dmoments = Moment[]
    for moment in reverse(moments)
        # backward step of quantum state
        (m, ρ, Δ) = backward_density(moment, ρ, Δ)
        pushfirst!(dmoments, m)
    end
    return dmoments, Δ
end


"""
    gradients(c::Circuit{N}, ρ::DensityMatrix, Δ::AbstractVector{<:Real}) where {N}

Perform a backward pass to compute gradients of a (fictitious) cost function with
respect to the circuit parameters of `c` and input density matrix `ρ`. `Δ` contains
the cost function gradients with respect to the circuit outputs (measurement operator averages).
The gradients with respect to the circuit parameters are returned in a duplicate circuit;
the overall return value is the tuple (dc::Circuit{N}, dρ::DensityMatrix).
"""
function gradients(c::Circuit{N}, ρ::DensityMatrix, Δ::AbstractVector{<:Real}) where {N}
    @assert(ρ.N == N)
    # length of circuit output vector must match gradient vector
    @assert(length(Δ) == length(c.meas))
    # forward pass through unitary gates
    apply!(ρ, c.moments)
    # gradient of cost function with respect to ρ
    ρbar = density_from_matrix(0.5^N * sum([Δ[i] * sparse_matrix(c.meas[i], N) for i in 1:length(Δ)]))
    # backward pass through unitary gates
    dmeas = MeasurementOperator.([Δ[i] * matrix(ρ) for i in 1:length(Δ)], (Tuple(1:N),))
    dcgc, ρbar = backward_density(c.moments, ρ, ρbar)
    return Circuit{N}(dcgc, dmeas), ρbar
end

"""
    gradients!(c::Circuit{N}, ρ::DensityMatrix, Δ::AbstractVector{<:Real}) where {N}

Perform a backward pass to compute gradients of a (fictitious) cost function with
respect to the circuit parameters of `c` and input density matrix `ρ`. `Δ` contains
the cost function gradients with respect to the circuit outputs (measurement operator averages).
The gradients with respect to the circuit parameters are returned in a duplicate circuit;
the overall return value is the tuple (dc::Circuit{N}, dρ::DensityMatrix). Custom helper function
works for backpropagation for in-place apply! method.
"""
function gradients!(c::Circuit{N}, ρ::DensityMatrix, Δ::AbstractVector{<:Real}) where {N}
    @assert(ρ.N == N)
    # length of circuit output vector must match gradient vector
    @assert(length(Δ) == length(c.meas))
    # gradient of cost function with respect to ρ
    ρbar = density_from_matrix(0.5^N * sum([Δ[i] * sparse_matrix(c.meas[i], N) for i in 1:length(Δ)]))
    # backward pass through unitary gates
    dmeas = MeasurementOperator.([Δ[i] * matrix(ρ) for i in 1:length(Δ)], (Tuple(1:N),))
    dcgc, ρbar = backward_density(c.moments, ρ, ρbar)
    return Circuit{N}(dcgc, dmeas), ρbar
end

