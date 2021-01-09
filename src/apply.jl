
"""
    apply(m::AbstractMatrix, ψ::AbstractVector)

"""
apply(m::AbstractMatrix, ψ::Vector{<:Complex}) = m*ψ

_apply(m::AbstractMatrix, ψ::Vector{<:Complex}) = m*ψ


"""
    apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}

Apply a `CircuitGate{M,G}` to a quantum state vector `ψ`.
"""
function apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")

    _apply(cg, ψ, N)
end

"""
    apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})

Apply a `CircuitGate{M,G}` to a quantum state vector `ψ`.
"""
function apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    ψl = deepcopy(ψ)
    l == 2^N || error("Vector length must be a power of 2")
    for cg in cgs 
        ψl = _apply(cg, ψl, N)
    end
    ψl
end

function _apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}, N::Int) where {M,G}
    gtuples = reshape(Qaintessent.cartesian_tuples(2, M), :)
    U = matrix(cg.gate)
    ψ = reshape(ψ, fill(2, N)...)

    ψs = similar(ψ)
    for (i, it) in enumerate(gtuples)
        # cannot use .= here since broadcasting fails for scalar numbers
        ψs[sliced_index(it, cg.iwire, N)...] = sum(U[i, j] .* ψ[sliced_index(jt, cg.iwire, N)...] for (j, jt) in enumerate(gtuples))
    end
    return reshape(ψs, :)
end



"""Tailored apply for XGate"""
function _apply(cg::CircuitGate{1,XGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    ψ = ψ[:, [2, 1], :]
    return reshape(ψ, :)
end


"""Tailored apply for YGate"""
function _apply(cg::CircuitGate{1,YGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)

    A[:, 1, :] .= -im.*ψ[:, 2, :]
    A[:, 2, :] .=  im.*ψ[:, 1, :]
    return reshape(A, :)
end


"""Tailored apply for ZGate"""
function _apply(cg::CircuitGate{1,ZGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .=   ψ[:, 1, :]
    A[:, 2, :] .= .-ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for HadamardGate"""
function _apply(cg::CircuitGate{1,HadamardGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= (ψ[:, 1, :] .+ ψ[:, 2, :]) ./ sqrt(2)
    A[:, 2, :] .= (ψ[:, 1, :] .- ψ[:, 2, :]) ./ sqrt(2)
    return reshape(A, :)
end


"""Tailored apply for SGate"""
function _apply(cg::CircuitGate{1,SGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= im .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for SdagGate"""
function _apply(cg::CircuitGate{1,SdagGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= -im .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for TGate"""
function _apply(cg::CircuitGate{1,TGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= Base.exp(im*π/4) .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for TdagGate"""
function _apply(cg::CircuitGate{1,TdagGate}, ψ::Vector{<:Complex}, N::Int)
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= Base.exp(-im*π/4) .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for RzGate"""
function _apply(cg::CircuitGate{1,RzGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    θ = cg.gate.θ[]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= Base.exp(-im*θ/2).*ψ[:, 1, :]
    A[:, 2, :] .= Base.exp( im*θ/2).*ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for PhaseShiftGate"""
function _apply(cg::CircuitGate{1,PhaseShiftGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    A[:, 1, :] .= ψ[:, 1, :]
    A[:, 2, :] .= Base.exp(im*cg.gate.ϕ[]) .* ψ[:, 2, :]
    return reshape(A, :)
end


"""Tailored apply for a general single qubit gate"""
function _apply(cg::CircuitGate{1,AbstractGate}, ψ::Vector{<:Complex}, N::Int)
    i = cg.iwire[1]
    ψ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    A = similar(ψ)
    U = matrix(cg.gate)
    A[:, 1, :] .= U[1, 1] .* ψ[:, 1, :] .+ U[1, 2] .* ψ[:, 2, :]
    A[:, 2, :] .= U[2, 1] .* ψ[:, 1, :] .+ U[2, 2] .* ψ[:, 2, :]
    return reshape(A, :)
end
#

"""Tailored apply for SwapGate"""
function _apply(cg::CircuitGate{2,SwapGate}, ψ::Vector{<:Complex}, N::Int)
    i,j = cg.iwire
    i,j = i<j ? (i,j) : (j,i) #sort them
    qubit_slices = [i-1, 1, j-i-1, 1, N-j]

    ψr = reshape(ψ, 2 .^qubit_slices...)
    A = similar(ψr)

    A[:, 1, :, 1, :] .= ψr[:, 1, :, 1, :]
    A[:, 2, :, 2, :] .= ψr[:, 2, :, 2, :]
    A[:, 1, :, 2, :] .= ψr[:, 2, :, 1, :]
    A[:, 2, :, 1, :] .= ψr[:, 1, :, 2, :]

    return reshape(A, :)
end


"""Tailored apply for a general ControlledGate"""
function _apply(cg::CircuitGate{M,ControlledGate{G}}, ψ::Vector{<:Complex}, N::Int) where {M,G}
    A = copy(ψ) # TODO: copy only what's needed
    A = reshape(A, fill(2, N)...)
    T = target_wires(cg.gate)
    itarget  = cg.iwire[1:T]
    icontrol = cg.iwire[T+1:M]

    new_iwire = Tuple(i - sum(icontrol .< i) for i in itarget)

    slice_target = Tuple(i in icontrol ? 2 : Colon() for i in 1:N)
    new_cg = CircuitGate(new_iwire, cg.gate.U)

    A[slice_target...] = reshape(apply(new_cg, reshape(A[slice_target...], :)),
                                 fill(2, N-M+T)...)
    return reshape(A, :)
end

"""
    apply(m::Moment, ψ::Vector{<:Complex}) 

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""
function apply(m::Moment, ψ::Vector{<:Complex}) 
    Nmoment = req_wires(m)
    l = length(ψ)::Int    
    Nψ = Qaintessent.intlog2(l)
    l == 2^Nψ || error("Vector length must be a power of 2")
    Nmoment <= Nψ || error("Moment affecting $Nmoment qubits applied to $Nψ qubits")
    _apply(m, ψ, Nψ)   
end

function apply(m::Vector{Moment}, ψ::Vector{<:Complex})
    Nmoment = maximum(req_wires.(m))
    l = length(ψ)::Int    
    Nψ = Qaintessent.intlog2(l)
    l == 2^Nψ || error("Vector length must be a power of 2")
    Nmoment <= Nψ || error("Moment affecting $Nmoment qubits applied to $Nψ qubits")
    for moment in m
        ψ = _apply(moment, ψ, Nψ)
    end
    return ψ
end

function _apply(m::Moment, ψ::Vector{<:Complex}, N::Int)
    for gate in m
        ψ = _apply(gate, ψ, N)
    end
    return ψ
end

"""
    apply(m::MeasurementOperator, ψ::Vector{<:Complex}) 

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""

function apply(m::MeasurementOperator{M,G}, ψ::Vector{<:Complex}) where {M,G<:AbstractGate}
    Nmoment = num_wires(m)
    l = length(ψ)::Int    
    Nψ = Qaintessent.intlog2(l)
    l == 2^Nψ || error("Vector length must be a power of 2")
    Nmoment <= Nψ || error("MeasurementOperator affecting $Nmoment qubits applied to $Nψ qubits")
    c = circuit_gate((m.iwire...), m.operator)
    _apply(c, ψ, Nψ)   
end

function apply(m::MeasurementOperator{M,G}, ψ::Vector{<:Complex}) where {M,G<:AbstractMatrix}
    Nmoment = num_wires(m)
    l = length(ψ)::Int    
    Nψ = Qaintessent.intlog2(l)
    l == 2^Nψ || error("Vector length must be a power of 2")
    Nmoment <= Nψ || error("MeasurementOperator affecting $Nmoment qubits applied to $Nψ qubits")
    apply(m.operator, ψ)   
end


"""
    apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}

returns list of expectation values from measurement operators in `c.meas` after applying circuit gates in `c.cgc` on state vector of `N` qubits `ψ`
"""
function apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}
    length(ψ)::Int == 2^N || error("Size of vector `ψ` must match Circuit size of $(2^N)")
    ψl = deepcopy(ψ)
    for moment in c.moments
        ψl = _apply(moment, ψl, N)
    end
    ψr = apply.(c.meas, (ψl,)) 
    return real.(dot.((ψl,), ψr))
end


(c::Circuit)(ψ) = apply(c, ψ)
