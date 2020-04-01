
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

# apply circuit gate to quantum state vector
function oldapply(cg::CircuitGate{M,N}, ψ::AbstractVector) where {M,N}
    # TODO: optimize (do not explicitly generate matrix)
    return matrix(cg) * ψ
end



function bitswap(n::Int, p1::Int, p2::Int, N::Int)
    if p1 == p2
        return n
    end
    p1̄ = N-p1
    p2̄ = N-p2
    bit1 = (n >> p1̄) & 1
    bit2 = (n >> p2̄) & 1
    x = (bit1 ⊻ bit2)
    x = (x << p1̄) | (x << p2̄)
    return n ⊻ x
end

ror(x::Int, k::Int) = (x >>> (0x3f & k)) | (x << (0x3f & -k))
rol(x::Int, k::Int) = (x << (0x3f & k)) | (x >>> (0x3f & -k))

function swap(w::AbstractArray, p1::Int, p2::Int)
    tmp = w[p1]
    w[p1] = w[p2]
    w[p2] = tmp
    return w
end

function apply(cg::CircuitGate{M,N}, ψ::AbstractVector) where {M,N}
    W = length(cg.iwire)
    wires = [i for i in cg.iwire]

    gmat = matrix(cg.gate)
    buffer = Array{Complex{Float64}, 1}(undef, 2^N)
    # Use indices starting from 0 for easier calculation
    indices = 0:2^(N)-1

    for i in 1:W
        target = N-W+i
        indices = bitswap.(indices, wires[i], target, N)
        if target in wires
            target_index = Base.findfirst(x -> x==target, wires)
            wires = swap(wires, i, target_index)
        end
    end

    indices = sortperm(indices)

    ψ = ψ[indices]

    # smaller O(S^2) for loop for repeated long Matrix-Vector operations
    # S is size of gmat
    S = size(gmat)[1]
    # Still not optimal, non-contiguous (?)
    for i in 1:S
        buffer[i:S:end] = sum(ψ[j:S:end] .* gmat[i,j] for j in 1:S)
    end

    indices = sortperm(indices)
    return buffer[indices]
end

Base.adjoint(cg::CircuitGate{M,N}) where {M,N} = CircuitGate{M,N}(cg.iwire, Base.adjoint(cg.gate))


single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate{1}, N::Integer) = CircuitGate{1,N}((iwire,), gate)


two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate{2}, N::Integer) = CircuitGate{2,N}((iwire1, iwire2), gate)


# single control and target wire
controlled_circuit_gate(icntrl::Integer, itarget::Integer, U::AbstractGate{1}, N::Integer) =
    controlled_circuit_gate((icntrl,), (itarget,), U, N)

# single control wire
controlled_circuit_gate(icntrl::Integer, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {M} =
    controlled_circuit_gate((icntrl,), itarget, U, N)

# single target wire
controlled_circuit_gate(icntrl::NTuple{K, <:Integer}, itarget::Integer, U::AbstractGate{1}, N::Integer)  where {K} =
    controlled_circuit_gate(icntrl, (itarget,), U, N)

function controlled_circuit_gate(icntrl::NTuple{K, <:Integer}, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {K,M}
    # consistency checks
    K + M ≤ N || error("Number of control and target wires must be smaller than overall number of wires.")
    length(intersect(icntrl, itarget)) == 0 || error("Control and target wires must be disjoint.")

    CircuitGate{K+M,N}((icntrl..., itarget...), ControlledGate{M,K+M}(U))
end


"""
    CircuitGateChain{N}

Chain of quantum circuit gates.
"""
mutable struct CircuitGateChain{N}
    gates::AbstractVector{<:AbstractCircuitGate{N}}
end

# unitary matrix representation of a sequence of circuit gates
function matrix(cgc::CircuitGateChain{N}) where {N}
    prod(Tuple(matrix(g) for g in reverse(cgc.gates)))
end

# apply circuit gate chain to quantum state vector
function apply(cgc::CircuitGateChain{N}, ψ::AbstractVector) where {N}
    for gate in cgc
        ψ = apply(gate, ψ)
    end
    return ψ
end

function oldapply(cgc::CircuitGateChain{N}, ψ::AbstractVector) where {M,N}
    # TODO: optimize (do not explicitly generate matrix)
    for gate in cgc
        ψ = oldapply(gate, ψ)
    end
    return ψ
end

Base.adjoint(cgc::CircuitGateChain{N}) where {N} = CircuitGateChain{N}(Base.adjoint.(reverse(cgc.gates)))

# make CircuitGateChain iterable and indexable
function Base.getindex(cgc::CircuitGateChain{N}, i::Integer) where {N}
    1 <= i <= Base.length(cgc.gates) || throw(BoundsError(S, i))
    return cgc.gates[i]
end

function Base.iterate(cgc::CircuitGateChain{N}, state=1) where {N}
    return state > Base.length(cgc.gates) ? nothing : (cgc[state], state+1)
end

# implement methods required for iteration
function Base.firstindex(cgc::CircuitGateChain{N}) where {N}
    return 1
end
function Base.lastindex(cgc::CircuitGateChain{N}) where {N}
    return Base.length(cgc.gates)
end


"""
    comm(A, B)

Matrix commutator [A, B].
"""
comm(A::AbstractMatrix, B::AbstractMatrix) = A*B - B*A


"""
    MeasurementOps{N}

Pairwise commuting measurement operators (Hermitian matrices).
"""
struct MeasurementOps{N}
    mops::AbstractVector{<:AbstractMatrix}

    function MeasurementOps{N}(mop::AbstractMatrix) where {N}
        # TODO: support general "qudits"
        d = 2
        # consistency checks
        size(mop) == (d^N, d^N) || error("Measurement operator must be a 2^N × 2^N matrix.")
        mop ≈ Base.adjoint(mop) || error("Measurement operator must be Hermitian.")
        new([mop])
    end

    function MeasurementOps{N}(mops::AbstractVector{<:AbstractMatrix}) where {N}
        # TODO: support general "qudits"
        d = 2
        # consistency checks
        for m in mops
            size(m) == (d^N, d^N) || error("Measurement operator must be a 2^N × 2^N matrix.")
            m ≈ Base.adjoint(m) || error("Measurement operator must be Hermitian.")
            for n in mops
                norm(comm(m, n))/d^N < 1e-13 || error("Measurement operators must pairwise commute.")
            end
        end
        new(mops)
    end
end


"""
    Circuit{N}

Quantum circuit consisting of a unitary gate chain and measurement operators.
"""
struct Circuit{N}
    cgc::CircuitGateChain{N}
    meas::MeasurementOps{N}
end


# apply circuit to quantum state vector and compute measurement expectation values
function apply(c::Circuit{N}, ψ::AbstractVector) where {N}
    ψs = apply(c.cgc, ψ)
    return [real(dot(ψs, m*ψs)) for m in c.meas.mops]
end

# get cost of quantum state vector and compute measurement expectation values
function cost(c::Circuit{N}, ψ::AbstractVector, e::Real) where {N}
    return sum([real(dot(ψ, m*ψ)) for m in c.meas.mops]) - e
end
