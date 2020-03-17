
"""
    AbstractCircuitGate{N}

Abtract unitary quantum circuit gate. `N` is the overall number of quantum "wires" of the circuit.
"""
abstract type AbstractCircuitGate{N} end


struct CircuitGate{M,N} <: AbstractCircuitGate{N}
    "ordered wire indices which this gate acts on"
    iwire::NTuple{M, <:Integer}
    "actual gate"
    gate::AbstractGate{M}

    function CircuitGate{M,N}(iwire::NTuple{M, <:Integer}, gate::AbstractGate{M}) where {M,N}
        M ≥ 1 || error("Need at least one wire to act on.")
        M ≤ N || error("Number of gate wires cannot be larger than total number of wires.")
        length(unique(iwire)) == M || error("Wire indices must be unique.")
        minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")
        maximum(iwire) ≤ N || error("Wire index cannot be larger than total number of wires.")
        new{M,N}(iwire, gate)
    end
end


function matrix(cg::CircuitGate{M,N}) where {M,N}
    # complementary wires
    iwcompl = setdiff(1:N, cg.iwire)
    @assert length(cg.iwire) + length(iwcompl) == N

    # TODO: support general "qudits"
    d = 2

    # TODO: handle sparse matrices efficiently
    gmat = matrix(cg.gate)
    @assert size(gmat) == (d^M, d^M)

    # Note: following the ordering convention of `kron` here, i.e.,
    # last qubit corresponds to fastest varying index
    strides = [d^(N-j) for j in 1:N]

    rowind = fill(0, d^(N+M))
    colind = fill(0, d^(N+M))
    values = fill(zero(eltype(gmat)), d^(N+M))
    count = 0
    for kw in (N-M > 0 ? reverse.(Iterators.product(fill(0:d-1, N-M)...)) : [()])
        for (i, iw) in enumerate(reverse.(Iterators.product(fill(0:d-1, M)...)))
            for (j, jw) in enumerate(reverse.(Iterators.product(fill(0:d-1, M)...)))
                # rearrange wire indices according to specification
                p = fill(0, N)
                for m in 1:(N-M); p[ iwcompl[m]] = kw[m]; end
                for m in 1:M;     p[cg.iwire[m]] = iw[m]; end
                q = fill(0, N)
                for m in 1:(N-M); q[ iwcompl[m]] = kw[m]; end
                for m in 1:M;     q[cg.iwire[m]] = jw[m]; end
                count += 1
                rowind[count] = dot(p, strides) + 1
                colind[count] = dot(q, strides) + 1
                values[count] = gmat[i, j]
            end
        end
    end
    @assert count == d^(N+M)

    return dropzeros!(sparse(rowind, colind, values, d^N, d^N))
end


single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate{1}, N::Integer) = CircuitGate{1,N}((iwire,), gate)


two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate{2}, N::Integer) = CircuitGate{2,N}((iwire1, iwire2), gate)


# single control and target wire
controlled_circuit_gate(icntrl::Integer, itarget::Integer, U::AbstractGate{1}, N::Integer) =
    controlled_circuit_gate((icntrl,), (itarget,), U, N)

# single control wire
controlled_circuit_gate(icntrl::Integer, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {M} =
    controlled_circuit_gate((icntrl,), itarget, U, N)

function controlled_circuit_gate(icntrl::NTuple{K, <:Integer}, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {K,M}
    # consistency checks
    K + M ≤ N || error("Number of control and target wires must be smaller than overall number of wires.")
    length(intersect(icntrl, itarget)) == 0 || error("Control and target wires must be disjoint.")

    CircuitGate{K+M,N}((icntrl..., itarget...), ControlledGate{M,K+M}(U))
end


"""
    CircuitBlock{N}

Quantum circuit block consisting of a sequence of gates.
"""
mutable struct CircuitBlock{N}
    gates::AbstractVector{<:AbstractCircuitGate{N}}
end

# convert circuit block to matrix
function matrix(b::CircuitBlock{N}) where {N}
    prod(Tuple(matrix(g) for g in reverse(b.gates)))
end

# apply circuit block to input vector of qubits
function apply(b::CircuitBlock{N}, ψ::AbstractVector) where {N}
    for gate in b
        ψ = matrix(gate) * ψ
    end
    return ψ
end

# make CircuitBlock iterable and indexable
function Base.getindex(b::CircuitBlock{N}, i::Int64) where {N}
    1 <= i <= Base.length(b.gates) || throw(BoundsError(S, i))
    return b.gates[i]
end

function Base.iterate(b::CircuitBlock{N}, state=1) where {N}
    return state > Base.length(b.gates) ? nothing : (b[state], state+1)
end

# implement methods required for iteration
function Base.firstindex(b::CircuitBlock{N}) where {N}
    return 1
end
function Base.lastindex(b::CircuitBlock{N}) where {N}
    return Base.length(b.gates)
end
