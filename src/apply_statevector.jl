using Base
using Base: invpermute!!

"""Modifies perm vector in Statevector object to flip qubit at position `wire`"""
function flipqubit!(ψ::Statevector, wire::Int)
    @simd for i in 1:length(ψ)
        @inbounds ψ.perm[i] = (ψ.perm[i] - 1) ⊻ (2^(wire-1)) + 1
    end
end

"""Permutes perm vector of length 2^N in Statevector object, shifting qubits in "wires" to the slowest running qubit for index n"""
function orderbit(n::Int, wires::NTuple{K,Int}, N::Int, control::NTuple{L,Int}) where {K,L}
    result = 0
    w = N-length(control)-length(wires)
    c = N-length(control)
    r = 0

    for i in 1:N
        bit = n >> (i-1) & 1
        if i in control
            result = result | (bit << c)
            c += 1
        elseif i in wires
            result = result | (bit << w)
            w += 1
        else
            result = result | (bit << r)
            r += 1
        end
    end
    return result
end

"""Permutes perm vector in Statevector object, shifting qubits in "wires" to the fastest running qubit"""
function orderqubit!(ψ::Statevector, wires::NTuple{K,Int}, control::NTuple{L,Int}=()) where {K,L}
    @simd for i in 1:length(ψ)
        @inbounds ψ.perm[orderbit(i-1, wires, ψ.N, control) + 1] = i
    end
end

"""Swaps bit at position `p1` and `p2` integer `n`"""
function swapbit(n::Int, p1::Int, p2::Int)
    # Move p1'th to rightmost side
    bit1 = (n >> p1) & 1
 
    # Move p2'th to rightmost side
    bit2 = (n >> p2) & 1
 
    # XOR the two bits
    x = (bit1 ⊻ bit2)
 
    # Put the xor bit back to their original positions
    x = (x << p1) | (x << p2)
 
    # XOR 'x' with the original number so that the
    # two sets are swapped
    result = n ⊻ x
    return result
end

"""Permutes perm vector in Statevector object, swapping qubit at `wire1` and `wire2`"""
function swapqubit!(ψ::Statevector, wire1::Int, wire2::Int)
    if wire1 == wire2
        return
    end
    @simd for i in 1:length(ψ)
        @inbounds ψ.perm[i] = swapbit(ψ.perm[i]-1, wire1-1, wire2-1) + 1
    end
end

function resetpermute!(ψ::Statevector)
    @inbounds ψ.perm .= 1:2^ψ.N
    return
end

"""Tailored apply for XGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,XGate})
    wire = cg.iwire[1]
    flipqubit!(ψ, wire)
    @inbounds ψ.vec .= getindex.(Ref(ψ.state),ψ.perm)
    @inbounds ψ.state .= ψ.vec
    resetpermute!(ψ)
    return
end


"""Tailored apply for YGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,YGate})
    wire = cg.iwire[1]
    bin = 2^(wire-1)

    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = -im*ψ[x]
        else
            ψ[x] = im*ψ[x]
        end
    end

    flipqubit!(ψ, wire)
    @inbounds ψ.vec .= getindex.(Ref(ψ.state),ψ.perm)
    @inbounds ψ.state .= ψ.vec
    resetpermute!(ψ)
    return
end


"""Tailored apply for ZGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,ZGate}) 
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    
    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = -ψ[x]
        end
    end
    return
end


"""Tailored apply for HadamardGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,HadamardGate}) 
    mid = 2^(ψ.N-1)
    ψ.state .*= convert(ComplexQ, 1/sqrt(2))

    orderqubit!(ψ, cg.iwire)
    @inbounds ψ.vec .= getindex.(Ref(ψ.state),ψ.perm)  
    @views begin
        @inbounds ψ.state .= ψ.vec
        @inbounds ψ.state[1:mid] .+= ψ.vec[mid+1:end]

        @inbounds ψ.state[mid+1:2mid] .*= -1
        @inbounds ψ.state[mid+1:2mid] .+= ψ.vec[1:mid]

        @inbounds invpermute!!(ψ.state, ψ.perm)
        resetpermute!(ψ)
    end
    return
end


"""Tailored apply for SGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,SGate}) 
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    
    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = im*ψ[x]
        end
    end
    return
end


"""Tailored apply for SdagGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,SdagGate}) 
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    
    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = -im*ψ[x]
        end
    end
    return
end


"""Tailored apply for TGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,TGate}) 
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    eπ4 = Base.exp(im*π/4)
    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = eπ4*ψ[x]
        end
    end
    return
end


"""Tailored apply for TdagGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,TdagGate})
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    eπ4 = Base.exp(-im*π/4)
    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = eπ4*ψ[x]
        end
    end
    return
end

"""Tailored apply for RzGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,RzGate}) 
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    θ = cg.gate.θ[]
    neπ2 = Base.exp(-im*θ/2)
    eπ2 = Base.exp(im*θ/2)
    for x in 1:2^ψ.N
        if (x-1) & bin != bin
            ψ[x] = neπ2*ψ[x]
        else
            ψ[x] = eπ2*ψ[x]
        end
    end
    return
end


"""Tailored apply for PhaseShiftGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,PhaseShiftGate}) 
    wire = cg.iwire[1]
    bin = 2^(wire-1)
    eϕ = Base.exp(im*cg.gate.ϕ[])
    for x in 1:2^ψ.N
        if (x-1) & bin == bin
            ψ[x] = eϕ*ψ[x]
        end
    end
    return
end


"""Tailored apply for a general single qubit gate"""
function _apply!(ψ::Statevector, cg::CircuitGate{1,T}) where {T<:AbstractGate}
    U = matrix(cg.gate)::Array{ComplexQ,2}
    orderqubit!(ψ, cg.iwire)
    @inbounds ψ.vec .= getindex.(Ref(ψ.state),ψ.perm)    
    mid = 2^(ψ.N-length(cg.iwire))
    @views begin
        @inbounds mul!(ψ.state[1:mid], U[1,1], ψ.vec[1:mid])
        @inbounds mul!(ψ.state[mid+1:end], U[2,1], ψ.vec[1:mid])
        @inbounds mul!(ψ.vec[1:mid], U[1,2], ψ.vec[mid+1:end])
        @inbounds ψ.state[1:mid] .+= ψ.vec[1:mid]
        @inbounds mul!(ψ.vec[1:mid], U[2,2], ψ.vec[mid+1:end])
        @inbounds ψ.state[mid+1:mid+mid] .+= ψ.vec[1:mid]
        @inbounds invpermute!!(ψ.state, ψ.perm)
        resetpermute!(ψ)
    end 
    return
end


"""Tailored apply for SwapGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{2,SwapGate})
    i, j = cg.iwire
    i, j = i < j ? (i, j) : (j, i) # sort them
    swapqubit!(ψ, i, j)
    @inbounds ψ.vec .= getindex.(Ref(ψ.state),ψ.perm)
    @inbounds ψ.state .= ψ.vec
    resetpermute!(ψ)
    return
end

"""Tailored apply for a general ControlledGate"""
function _apply!(ψ::Statevector, cg::CircuitGate{M,ControlledGate{G}}) where {M,G}
    T = target_wires(cg.gate)
    C = control_wires(cg.gate)
    orderqubit!(ψ, cg.iwire[1:T], cg.iwire[T+1:end])
    @inbounds ψ.vec .= getindex.(Ref(ψ.state), ψ.perm)
    new_cg = CircuitGate((ψ.N-C:-1:ψ.N-C-T+1...,), cg.gate.U)

    start = length(ψ) - 2^(ψ.N-C) + 1
    ψ.vec[start:end] .= _apply(ψ.vec[start:end], new_cg, ψ.N-C)
    @inbounds ψ.state .= ψ.vec
    @inbounds invpermute!!(ψ.state, ψ.perm)
    resetpermute!(ψ)
    return
end


"""
    apply(ψ::Statevector, cg::CircuitGate{M,G}) where {M,G}

Apply a [`CircuitGate`](@ref) to a quantum state vector `ψ`.
# Examples

```jldoctest
julia> cg = circuit_gate(1, HadamardGate());
julia> ψ = Statevector(ComplexQ[1 0]);
julia> apply!(ψ, cg)
julia> ψ.state
2-element Array{Complex{Float64},1}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
function apply!(ψ::Statevector, cg::CircuitGate{M,G}) where {M,G}
    req_wires(cg) <= ψ.N || error("CircuitGate requires a minimum of $(req_wires(cg)) qubits, input vector `ψ` has $ψ.N qubits")
    _apply!(ψ, cg)
end

function _apply!(ψ::Statevector, cg::CircuitGate{M,G}) where {M,G}
    U = matrix(cg.gate)::Array{ComplexQ,2}
    
    T = target_wires(cg.gate)
    orderqubit!(ψ, cg.iwire[1:T], cg.iwire[T+1:end])
    
    @inbounds ψ.vec .= getindex.(Ref(ψ.state),ψ.perm)    
    step = 2^(ψ.N-length(cg.iwire))
    @views begin
        for i in 1:2^M
            @inbounds mul!(ψ.state[(i-1)*step+1:i*step], U[i,1], ψ.vec[1:step])
        end
        for j in 2:2^M
            for i in 1:2^M
                @inbounds mul!(ψ.vec[1:step], U[i,j], ψ.vec[(j-1)*step+1:j*step])
                @inbounds ψ.state[(i-1)*step+1:i*step] .+= ψ.vec[1:step]
            end
        end
        @inbounds invpermute!!(ψ.state, ψ.perm)
        resetpermute!(ψ)
    end
end


"""
    apply(ψ::Statevector, cgs::Vector{<:CircuitGate})

Apply a sequence of [`CircuitGate`](@ref)(s) to a quantum state vector `ψ`.

# Examples

```jldoctest
julia> cgs = [circuit_gate(1, HadamardGate()),
                circuit_gate(1, X),
                circuit_gate(1, Y)];
julia> ψ = Statevector(ComplexQ[1 0]);
julia> apply!(ψ, cgs)
julia> ψ.state
2-element Array{Complex{Float64},1}:
 0.0 - 0.7071067811865475im
 0.0 + 0.7071067811865475im
```
"""
function apply!(ψ::Statevector, cgs::Vector{<:CircuitGate})
    maximum(req_wires.(cgs)) <= ψ.N || error("CircuitGates require a minimum of $req qubits, input vector `ψ` has $N qubits")
    for cg in cgs 
        _apply!(ψ, cg)
    end
end

"""
    apply(ψ::Statevector, m::Moment)

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""
function apply!(ψ::Statevector, m::Moment) 
    Nmoment = req_wires(m)
    Nmoment <= ψ.N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    _apply!(ψ, m)   
end

function _apply!(ψ::Statevector, m::Moment)
    for gate in m
        _apply!(ψ, gate)
    end
end

function apply!(ψ::Statevector, m::Vector{Moment})
    length(m) != 0 || error("Vector of length 0 cannot be applied")
    Nmoment = maximum(req_wires.(m))
    Nmoment <= ψ.N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    for moment in m
        _apply!(ψ, moment)
    end
end

"""
    apply(ψ::Statevector, c::Circuit{N}) where {N}

returns list of expectation values from measurement operators in `c.meas` after applying circuit gates in `c.cgc` on state vector of `N` qubits `ψ`
"""
function apply!(ψ::Statevector, c::Circuit{N}) where {N}
    N == ψ.N || error("Size of vector `ψ` must match Circuit size of $(2^N)")
    length(c) != 0 || error("Circuit does not contain any gates")
    length(c.meas) != 0 || error("Circuit does not contain any measurement operators")
    for moment in c.moments
        _apply!(ψ, moment)
    end
    ψr = apply.((ψ.state,), c.meas)
    return real.(dot.((ψ.state,), ψr))
end
