
apply(ψ::Vector{<:Complex}, m::AbstractMatrix) = m*ψ

"""Tailored apply for XGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,XGate}, N::Int) 
    i = cg.iwire[1]
    χ = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = χ[:, [2, 1], :]
    return reshape(χ, :)
end


"""Tailored apply for YGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,YGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= -im.*ψs[:, 2, :]
    χ[:, 2, :] .=  im.*ψs[:, 1, :]
    return reshape(χ, :)
end


"""Tailored apply for ZGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,ZGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .=   ψs[:, 1, :]
    χ[:, 2, :] .= .-ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for HadamardGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,HadamardGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= (ψs[:, 1, :] .+ ψs[:, 2, :]) ./ sqrt(2)
    χ[:, 2, :] .= (ψs[:, 1, :] .- ψs[:, 2, :]) ./ sqrt(2)
    return reshape(χ, :)
end


"""Tailored apply for SGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,SGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= im .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for SdagGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,SdagGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= -im .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for TGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,TGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp(im*π/4) .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for TdagGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,TdagGate}, N::Int)
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp(-im*π/4) .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for RzGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,RzGate}, N::Int) 
    i = cg.iwire[1]
    θ = cg.gate.θ[]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= Base.exp(-im*θ/2).*ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp( im*θ/2).*ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for PhaseShiftGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,PhaseShiftGate}, N::Int) 
    i = cg.iwire[1]
    ψs = reshape(ψ, 2^(i-1), 2, 2^(N-i))
    χ = similar(ψs)
    χ[:, 1, :] .= ψs[:, 1, :]
    χ[:, 2, :] .= Base.exp(im*cg.gate.ϕ[]) .* ψs[:, 2, :]
    return reshape(χ, :)
end


"""Tailored apply for a general single qubit gate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{1,T}, N::Int) where {T<:AbstractGate}
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
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{2,SwapGate}, N::Int)
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
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{M,ControlledGate{G}}, N::Int) where {M,G}
    χ = copy(ψ)::Vector{ComplexQ} # TODO: copy only what's needed
    χ = reshape(χ, fill(2, N)...)
    T = target_wires(cg.gate)
    itarget  = cg.iwire[1:T]
    icontrol = cg.iwire[T+1:M]

    new_iwire = Tuple(i - sum(icontrol .< i) for i in itarget)

    slice_target = Tuple(i in icontrol ? 2 : Colon() for i in 1:N)
    new_cg = CircuitGate(new_iwire, cg.gate.U)

    χ[slice_target...] = reshape(apply(reshape(χ[slice_target...], :), new_cg),
                                 fill(2, N-M+T)...)
    return reshape(χ, :)
end


"""Tailored apply for a general CircuitGate"""
function _apply(ψ::Vector{<:Complex}, cg::CircuitGate{M,G}, N::Int) where {M,G}
    U = matrix(cg.gate)
    ψs = reshape(ψ, fill(2, N)...)
    χ = similar(ψs)
    it = binary_digits(M, 0)
    m = binary_digits(M, 0)
    for i in 1:2^M
        binary_digits!(it, i - 1)
        # cannot use .= here since broadcasting fails for scalar numbers
        χ[sliced_index(it, cg.iwire, N)...] = sum(U[i, j] .* ψs[sliced_index(binary_digits!(m, j - 1), cg.iwire, N)...] for j in 1:2^M)
    end
    return reshape(χ, :)
end


"""
    apply(ψ::Vector{<:Complex}, cg::CircuitGate{M,G}) where {M,G}

Apply a [`CircuitGate`](@ref) to a quantum state vector `ψ`.
# Examples

```jldoctest
julia> cg = circuit_gate(1, HadamardGate());
julia> ψ = [1 0];
julia> apply(ψ, cg)
2-element Array{Complex{FloatQ},1}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
function apply(ψ::Vector{<:Complex}, cg::CircuitGate{M,G}) where {M,G}
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    req_wires(cg) <= N || error("CircuitGate requires a minimum of $(req_wires(cg)) qubits, input vector `ψ` has $N qubits")
    _apply(ψ, cg, N)
end



"""
    apply(ψ::Vector{<:Complex}, cgs::Vector{<:CircuitGate})

Apply a sequence of [`CircuitGate`](@ref)(s) to a quantum state vector `ψ`.

# Examples

```jldoctest
julia> cgs = [circuit_gate(1, HadamardGate()),
                circuit_gate(1, X),
                circuit_gate(1, Y)];
julia> ψ = [1 0];
julia> apply(ψ, cgs)
2-element Array{Complex{FloatQ},1}:
 0.0 - 0.7071067811865475im
 0.0 + 0.7071067811865475im
```
"""
function apply(ψ::Vector{<:Complex}, cgs::Vector{<:CircuitGate})
    length(cgs) != 0 || error("Vector of length 0 cannot be applied")
    l = length(ψ)
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    req = maximum(req_wires.(cgs))
    req <= N || error("CircuitGates require a minimum of $req qubits, input vector `ψ` has $N qubits")
    for cg in cgs 
        ψ = _apply(ψ, cg, N)
    end
    return ψ
end


"""
    apply(ψ::Vector{<:Complex}, m::Moment)

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""
function apply(ψ::Vector{<:Complex}, m::Moment) 
    Nmoment = req_wires(m)
    l = length(ψ)
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    _apply(ψ, m, N)   
end

function _apply(ψ::Vector{<:Complex}, m::Moment, N::Int)
    for gate in m
        ψ = _apply(ψ, gate, N)
    end
    return ψ
end


function apply(ψ::Vector{<:Complex}, m::Vector{Moment})
    length(m) != 0 || error("Vector of length 0 cannot be applied")
    Nmoment = maximum(req_wires.(m))
    l = length(ψ)
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    for moment in m
        ψ = _apply(ψ, moment, N)
    end
    return ψ
end


"""
    apply(ψ::Vector{<:Complex}, m::MeasurementOperator{M,G}) where {M,G<:AbstractGate}

returns state vector of `N` qubits after applying a `MeasurementOperator` object to a quantum state vector of `N` qubits `ψ`
"""
function apply(ψ::Vector{<:Complex}, m::MeasurementOperator{M,G}) where {M,G<:AbstractGate}
    Nmoment = num_wires(m)
    l = length(ψ)::Int    
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    Nmoment <= N || error("MeasurementOperator affecting $Nmoment qubits applied to $N qubits")
    c = circuit_gate((m.iwire...), m.operator)
    _apply(ψ, c, N)   
end


function apply(ψ::Vector{<:Complex}, m::MeasurementOperator{M,G}) where {M,G<:AbstractMatrix}
    l = length(ψ)::Int
    N = Qaintessent.intlog2(l)
    l == 2^N || error("Vector length must be a power of 2")
    M <= N || error("MeasurementOperator affecting $M qubits applied to $N qubits")
    U = m.operator
    ψs = reshape(ψ, fill(2, N)...)
    χ = similar(ψs)
    it = binary_digits(M, 0)
    scratch = binary_digits(M, 0)
    for i in 1:2^M
        binary_digits!(it, i - 1)
        # cannot use .= here since broadcasting fails for scalar numbers
        χ[sliced_index(it, m.iwire, N)...] = sum(U[i, j] .* ψs[sliced_index(binary_digits!(scratch, j - 1), m.iwire, N)...] for j in 1:2^M)
    end
    return reshape(χ, :)
end


"""
    apply(ψ::Vector{<:Complex}, c::Circuit{N}) where {N}

returns list of expectation values from measurement operators in `c.meas` after applying circuit gates in `c.cgc` on state vector of `N` qubits `ψ`
"""
function apply(ψ::Vector{<:Complex}, c::Circuit{N}) where {N}
    length(ψ)::Int == 2^N || error("Size of vector `ψ` must match Circuit size of $(2^N)")
    length(c) != 0 || error("Circuit does not contain any gates")
    length(c.meas) != 0 || error("Circuit does not contain any measurement operators")
    ψl = copy(ψ)
    for moment in c.moments
        ψl = _apply(ψl, moment, N)
    end
    ψr = apply.((ψl,), c.meas) 
    return real.(dot.((ψl,), ψr))
end

(c::Circuit)(ψ) = apply(ψ, c)