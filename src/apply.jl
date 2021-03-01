
apply(m::AbstractMatrix, ψ::Vector{<:Complex}) = m*ψ

_apply(m::AbstractMatrix, ψ::Vector{<:Complex}) = m*ψ



"""Tailored apply for XGate"""
function _apply(cg::CircuitGate{1,XGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    χ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = χ[:, [2, 1], :]
    return reshape(χ, :)
end


"""Tailored apply for YGate"""
function _apply(cg::CircuitGate{1,YGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= -im.*ψs[:, 2, :]
    χ[:, 2, :] .=  im.*ψs[:, 1, :]
    return reshape(χ, :)
end


"""Tailored apply for ZGate"""
function _apply(cg::CircuitGate{1,ZGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .=   ψs[:, 1, :]
    χ[:, 2, :] .= .-ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for HadamardGate"""
function _apply(cg::CircuitGate{1,HadamardGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= (ψs[:, 1, :] .+ ψs[:, 2, :]) ./ sqrt(2)
    χ[:, 2, :] .= (ψs[:, 1, :] .- ψs[:, 2, :]) ./ sqrt(2)
    return reshape(χ, :)
end


"""Tailored apply for SGate"""
function _apply(cg::CircuitGate{1,SGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= im .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for SdagGate"""
function _apply(cg::CircuitGate{1,SdagGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= -im .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for TGate"""
function _apply(cg::CircuitGate{1,TGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp(im*π/4) .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for TdagGate"""
function _apply(cg::CircuitGate{1,TdagGate}, ψ::Vector{<:Complex}, N::Int)
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp(-im*π/4) .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for RzGate"""
function _apply(cg::CircuitGate{1,RzGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    θ = cg.gate.θ[]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= Base.exp(-im*θ/2).*ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp( im*θ/2).*ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for PhaseShiftGate"""
function _apply(cg::CircuitGate{1,PhaseShiftGate}, ψ::Vector{<:Complex}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp(im*cg.gate.ϕ[]) .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for a general single qubit gate"""
function _apply(cg::CircuitGate{1,AbstractGate}, ψ::Vector{<:Complex}, N::Int)
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    U = matrix(cg.gate)
    χ[:, 1, :] .= U[1, 1] .* ψs[:, 1, :] .+ U[1, 2] .* ψs[:, 2, :]
    χ[:, 2, :] .= U[2, 1] .* ψs[:, 1, :] .+ U[2, 2] .* ψs[:, 2, :]
    return reshape(χ, :)
end
#

"""Tailored apply for SwapGate"""
function _apply(cg::CircuitGate{2,SwapGate}, ψ::Vector{<:Complex}, N::Int)
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i) # sort them

    ψs = reshape(ψ, 2^(i-1), 2, 2^(j-i-1), 2, 2^(N-j))
    χ = similar(ψs)

    χ[:, 1, :, 1, :] .= ψs[:, 1, :, 1, :]
    χ[:, 2, :, 2, :] .= ψs[:, 2, :, 2, :]
    χ[:, 1, :, 2, :] .= ψs[:, 2, :, 1, :]
    χ[:, 2, :, 1, :] .= ψs[:, 1, :, 2, :]

    return reshape(χ, :)
end


"""Tailored apply for a general ControlledGate"""
function _apply(cg::CircuitGate{M,ControlledGate{G}}, ψ::Vector{<:Complex}, N::Int) where {M,G}
    χ = copy(ψ) # TODO: copy only what's needed
    χ = reshape(χ, fill(2, N)...)
    T = target_wires(cg.gate)
    itarget  = cg.iwire[1:T]
    icontrol = cg.iwire[T+1:M]

    new_iwire = Tuple(i - sum(icontrol .< i) for i in itarget)

    slice_target = Tuple(i in icontrol ? 2 : Colon() for i in 1:N)
    new_cg = CircuitGate(new_iwire, cg.gate.U)

    χ[slice_target...] = reshape(apply(new_cg, reshape(χ[slice_target...], :)),
                                 fill(2, N-M+T)...)
    return reshape(χ, :)
end


"""
    apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}

Apply a [`CircuitGate`](@ref) to a quantum state vector `ψ`.
# Examples

```jldoctest
julia> cg = circuit_gate(1, HadamardGate());
julia> ψ = [1 0];
julia> apply(cg, ψ)
2-element Array{Complex{Float64},1}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
function apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}) where {M,G}
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    req_wires(cg) <= N || error("CircuitGate requires a minimum of $(req_wires(cg)) qubits, input vector `ψ` has $N qubits")
    _apply(cg, ψ, N)
end

function _apply(cg::CircuitGate{M,G}, ψ::Vector{<:Complex}, N::Int) where {M,G}
    U = matrix(cg.gate)
    ψs = reshape(ψ, fill(2, N)...)
    χ = similar(ψs)
    for i in 1:2^M
        it = binary_digits(i - 1, M)
        # cannot use .= here since broadcasting fails for scalar numbers
        χ[sliced_index(it, cg.iwire, N)...] = sum(U[i, j] .* ψs[sliced_index(binary_digits(j - 1, M), cg.iwire, N)...] for j in 1:2^M)
    end
    return reshape(χ, :)
end


"""
    apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})

Apply a sequence of [`CircuitGate`](@ref)(s) to a quantum state vector `ψ`.

# Examples

```jldoctest
julia> cgs = [circuit_gate(1, HadamardGate()),
                circuit_gate(1, X),
                circuit_gate(1, Y)];
julia> ψ = [1 0];
julia> apply(cgs, ψ)
2-element Array{Complex{Float64},1}:
 0.0 - 0.7071067811865475im
 0.0 + 0.7071067811865475im
```
"""
function apply(cgs::Vector{<:CircuitGate}, ψ::Vector{<:Complex})
    l = length(ψ)
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    req = maximum(req_wires.(cgs))
    req <= N || error("CircuitGates require a minimum of $req qubits, input vector `ψ` has $N qubits")
    for cg in cgs 
        ψ = _apply(cg, ψ, N)
    end
    return ψ
end


"""
    apply(m::Moment, ψ::Vector{<:Complex}) 

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""
function apply(m::Moment, ψ::Vector{<:Complex}) 
    Nmoment = req_wires(m)
    l = length(ψ)
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    _apply(m, ψ, N)   
end

function _apply(m::Moment, ψ::Vector{<:Complex}, N::Int)
    for gate in m
        ψ = _apply(gate, ψ, N)
    end
    return ψ
end


function apply(m::Vector{Moment}, ψ::Vector{<:Complex})
    Nmoment = maximum(req_wires.(m))
    l = length(ψ)
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    for moment in m
        ψ = _apply(moment, ψ, N)
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
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("MeasurementOperator affecting $Nmoment qubits applied to $N qubits")
    c = circuit_gate((m.iwire...), m.operator)
    _apply(c, ψ, N)   
end

function apply(m::MeasurementOperator{M,G}, ψ::Vector{<:Complex}) where {M,G<:AbstractMatrix}
    Nmoment = num_wires(m)
    l = length(ψ)::Int    
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("MeasurementOperator affecting $Nmoment qubits applied to $N qubits")
    apply(m.operator, ψ)   
end


"""
    apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}

returns list of expectation values from measurement operators in `c.meas` after applying circuit gates in `c.cgc` on state vector of `N` qubits `ψ`
"""
function apply(c::Circuit{N}, ψ::Vector{<:Complex}) where {N}
    length(ψ)::Int == 2^N || error("Size of vector `ψ` must match Circuit size of $(2^N)")
    ψl = copy(ψ)
    for moment in c.moments
        ψl = _apply(moment, ψl, N)
    end
    ψr = apply.(c.meas, (ψl,)) 
    return real.(dot.((ψl,), ψr))
end


(c::Circuit)(ψ) = apply(c, ψ)
