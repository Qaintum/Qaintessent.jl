
"""
    AbstractCircuitGate{N}

Abstract unitary quantum circuit gate. `N` is the overall number of quantum "wires" of the circuit.
"""
abstract type AbstractCircuitGate{N} end


function int2bit(x::Int64; pad=nothing)
    if isnothing(pad)
        pad = floor(Int64, log2(3)+1)
    end
    BitArray(digits(x, base=2, pad=pad) |> reverse)
end

function bit2int(b::BitArray)
    b = deepcopy(b)
    num = 0
    count = 0
    while !isempty(b)
        num += pop!(b)*2^count
        count += 1
    end
    num
end

"""
    CircuitGate{M,N,G} <: AbstractCircuitGate{N}

Unitary quantum circuit gate. `M` is the number of wires affected by the CircuitGate, `N` is the overall number of quantum "wires" of the circuit, `G` is the basic gate used to construct the CircuitGate.
"""
struct CircuitGate{M,N,G} <: AbstractCircuitGate{N}
    "ordered wire indices which this gate acts on"
    iwire::NTuple{M,<:Integer}
    "actual gate"
    gate::G
    "classical registers"
    ccntrl::AbstractVector{<:Union{Int, Expr}}

    @doc """
        CircuitGate{M,N,G}(iwire::NTuple{M,<:Integer}, gate::G) where {M,N,G}

    Creates a `CircuitGate{M,N,G}` object. `M` is the number of wires affected by the CircuitGate, `N` is the overall number of quantum "wires" of the circuit, `G` is the basic gate used to construct the CircuitGate.
    """
    function CircuitGate{M,N,G}(iwire::NTuple{M, <:Integer}, gate::G; ccntrl::AbstractVector{<:Union{Int, Expr}}=Int[]) where {M,N,G}
        M ≥ 1 || error("Need at least one wire to act on.")
        M ≤ N || error("Number of gate wires cannot be larger than total number of wires.")
        length(unique(iwire)) == M || error("Wire indices must be unique.")
        minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")
        maximum(iwire) ≤ N || error("Wire index cannot be larger than total number of wires.")
        G <: AbstractGate{M} || error("Gate type must be a subtype of AbstractGate{$M}.")

        if ccntrl isa AbstractVector{Integer} && !isempty(ccntrl)
            sum(ccntrl .> 0) == length(ccntrl) || error("Classical registers cannot be negative and indexing starts from 1")
        end

        new{M,N,G}(iwire, gate, ccntrl)
    end
end

"""
    CircuitGate(iwire::NTuple{M,<:Integer}, gate::AbstractGate{M}, N) where {M}

Creates a `CircuitGate{M,N,G}` object. `M` is the number of wires affected by the CircuitGate, `N` is the overall number of quantum "wires" of the circuit, `G` is the basic gate used to construct the CircuitGate.
"""
function CircuitGate(iwire::NTuple{M, <:Integer}, gate::AbstractGate{M}, N; ccntrl::Array{<:Union{Int64,Expr},1}=Int[]) where {M}
    G = typeof(gate)
    CircuitGate{M,N,G}(iwire, gate; ccntrl=ccntrl)
end

"""
    Base.isapprox(cg1::CircuitGate{M,N,G}, cg2::CircuitGate{M,N,G})

Compares two circuit gates of basic type `G`.
"""
function Base.isapprox(cg1::CircuitGate{M,N,G}, cg2::CircuitGate{M,N,G}) where {M,N,G}
    # for parametric gates, return true only if parameters are approximately equal
    fields = fieldnames(G)
    if cg1.iwire != cg2.iwire
        return false
    end
    for name in fields
        if !(getfield(cg1.gate, name) ≈ getfield(cg2.gate, name))
            return false
        end
    end
    return true
end

Base.isapprox(cg1::CircuitGate{M, N, G}, cg2::CircuitGate{M, N, H}) where {M, N, G, H} = false

LinearAlgebra.ishermitian(cg::CircuitGate) = LinearAlgebra.ishermitian(cg.gate)

"""
    matrix(cg::CircuitGate{M,N,G}) where {M,N,G<:AbstractGate}

returns matrix representation of a `CircuitGate{M,N,G}` object that can applied to a state vector of `N` qubits.
"""
function matrix(cg::CircuitGate{M,N,G}) where {M,N,G<:AbstractGate}
    # convert to array
    iwire = collect(cg.iwire)
    # complementary wires
    iwcompl = setdiff(1:N, iwire)
    @assert length(iwire) + length(iwcompl) == N

    # TODO: support general "qudits"
    d = 2

    # TODO: handle sparse matrices efficiently
    gmat = matrix(cg.gate)
    @assert size(gmat) == (d^M, d^M)

    # Note: following the ordering convention of `kron` here, i.e.,
    # last qubit corresponds to fastest varying index
    strides = [d^(j-1) for j in 1:N]
    wstrides = strides[iwire]
    cstrides = strides[iwcompl]

    rowind = fill(0, d^(N+M))
    colind = fill(0, d^(N+M))
    values = fill(zero(eltype(gmat)), d^(N+M))
    count = 0
    for kw in (N-M > 0 ? reverse.(cartesian_tuples(d, N-M)) : [Int[]])
        koffset = dot(collect(kw), cstrides)
        for (i, iw) in enumerate(reverse.(cartesian_tuples(d, M)))
            for (j, jw) in enumerate(reverse.(cartesian_tuples(d, M)))
                # rearrange wire indices according to specification
                count += 1
                rowind[count] = koffset + dot(collect(iw), wstrides) + 1
                colind[count] = koffset + dot(collect(jw), wstrides) + 1
                values[count] = gmat[i, j]
            end
        end
    end
    @assert count == d^(N+M)

    return dropzeros!(sparse(rowind, colind, values, d^N, d^N))
end

"""
    Base.adjoint(cg::CircuitGate{M,N,G})

Construct a `CircuitGate{M,N,H}` object where `H` is the adjoint of `AbstractGate` `G`
"""
function Base.adjoint(cg::CircuitGate{M,N,G}) where {M,N,G}
    adj_gate = Base.adjoint(cg.gate)
    CircuitGate{M,N,typeof(adj_gate)}(cg.iwire, adj_gate)
end

"""
    single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate{1}, N::Integer)

Construct a `CircuitGate{1,N,G}` object of basic gate type `gate` affecting wire `iwire`.
"""
single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate{1}, N::Integer) =
    CircuitGate((iwire,), gate, N)

"""
    single_qubit_circuit_gate(qreg::QRegister, gate::AbstractGate{1}, N::Integer)

Construct a `CircuitGate{1,N,G}` objects of basic gate type `gate` affecting QRegister `qreg`
"""
function single_qubit_circuit_gate(qreg::QRegister, gate::AbstractGate{1}, N::Integer)
    !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    [single_qubit_circuit_gate(i, gate, N) for i in qreg.ind]
end

"""
    two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate{2}, N::Integer)

Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting wires `iwire1` and `iwire2`.
"""
two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate{2}, N::Integer) =
    CircuitGate((iwire1, iwire2), gate, N)

"""
    two_qubit_circuit_gate(qreg::QRegister, iwire2::Integer, gate::AbstractGate{2}, N::Integer)

Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting QRegister `qreg` and wire `iwire2`.
"""
function two_qubit_circuit_gate(qreg::QRegister, iwire2::Integer, gate::AbstractGate{2}, N::Integer)
    !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    [two_qubit_circuit_gate(i, iwire2, gate, N) for i in qreg.ind]
end

"""
    two_qubit_circuit_gate(iwire1::Integer, qreg::QRegister, gate::AbstractGate{2}, N::Integer)

Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting QRegister `qreg` and wire `iwire1`.
"""
function two_qubit_circuit_gate(iwire1::Integer, qreg::QRegister, gate::AbstractGate{2}, N::Integer)
    !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    [two_qubit_circuit_gate(iwire1, i, gate, N) for i in qreg.ind]
end

"""
    two_qubit_circuit_gate(qreg1::QRegister, qreg1::QRegister, gate::AbstractGate{2}, N::Integer)

Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting QRegister `qreg1` and ``qreg2`.
"""
function two_qubit_circuit_gate(qreg1::QRegister, qreg2::QRegister, gate::AbstractGate{2}, N::Integer)
    !any(qreg1.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    !any(qreg2.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    length(qreg1.ind) == length(qreg2.ind) || error("Only able to apply CircuitGate to two QRegister objects if registers are of same length")
    [two_qubit_circuit_gate(i, j, gate, N) for (i,j) in zip(qreg1.ind, qreg2.ind)]
end

# single control and target wire
"""
    controlled_circuit_gate(icntrl::Union{Integer, Expr}, itarget::Integer, U::AbstractGate{1}, N::Integer)

Construct a `CircuitGate{2,N,G}` object of basic gate type `U` controlled by wire or Expr `icntrl` and affecting wire `itarget`.
"""
controlled_circuit_gate(icntrl::Union{Integer, Expr}, itarget::Integer, U::AbstractGate{1}, N::Integer) =
    controlled_circuit_gate((icntrl,), (itarget,), U, N)


# single control wire
"""
    controlled_circuit_gate(icntrl::Integer, itarget::NTuple{M,<:Integer}, U::AbstractGate{M}, N::Integer) where {M}

Construct a `CircuitGate{M+1,N,G}` object of basic gate type `U` controlled by wire or Expr `icntrl` and affecting wires in tuple `itarget`.
"""
controlled_circuit_gate(icntrl::Union{Integer, Expr}, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {M} =
    controlled_circuit_gate((icntrl,), itarget, U, N)

# single target wire
"""
    controlled_circuit_gate(icntrl::NTuple{K, Union{Int, Expr}}, itarget::Integer, U::AbstractGate{1}, N::Integer)  where {K}

Construct a `CircuitGate{K+1,N,G}` object of basic gate type `U` controlled by wires or Expr in tuple `icntrl` and affecting wire `itarget`.
"""
controlled_circuit_gate(icntrl::NTuple{K, Union{Integer, Expr}}, itarget::Integer, U::AbstractGate{1}, N::Integer)  where {K} =
    controlled_circuit_gate(icntrl, (itarget,), U, N)

"""
    controlled_circuit_gate(icntrl::NTuple{K, <:Union{Integer, Expr}}, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {K,M}

Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in tuple `icntrl` and affecting wires in tuple `itarget`.
"""
function controlled_circuit_gate(icntrl::NTuple{K, Union{Integer, Expr}}, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer) where {K,M}
    ccntrl = Union{Int, Expr}[]
    for wire in icntrl
        if wire isa Expr
            push!(ccntrl, wire)
        elseif wire < 0
            push!(ccntrl, abs(wire))
        end
    end
    icntrl = (Iterators.filter(x->x isa Int && x>0, icntrl)...,)
    length(ccntrl) == length(unique(ccntrl)) || error("Classical control wires must be unique")

    if length(icntrl) == 0
        return CircuitGate((itarget...,), U, N; ccntrl=ccntrl)
    end

    # consistency checks
    k = length(icntrl)
    k + M ≤ N || error("Number of control and target wires must be smaller than overall number of wires.")
    length(intersect(icntrl, itarget)) == 0 || error("Control and target wires must be disjoint.")

    CircuitGate((icntrl..., itarget...), ControlledGate{M,k+M}(U), N; ccntrl=ccntrl)
end

# control with QRegisters
"""
    controlled_circuit_gate(reg::Register, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer; ccntrl::AbstractVector{Integer}=Int[]) where {K,M}

Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in register `reg` and affecting wires in tuple `itarget`.
"""
function controlled_circuit_gate(reg::Register, itarget::NTuple{M, <:Integer}, U::AbstractGate{M}, N::Integer; ccntrl::AbstractVector{<:Integer}=Int[]) where {K,M}
    !any(reg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    [controlled_circuit_gate(i, itarget, U, N) for i in reg.ind]
end

controlled_circuit_gate(reg::Register, itarget::Integer, U::AbstractGate{M}, N::Integer; ccntrl::AbstractVector{<:Integer}=Int[]) where {K,M} =
    controlled_circuit_gate(reg, (itarget,), U, N)


"""
    controlled_circuit_gate(qreg::QRegister, itarget::NTuple{M,<:Integer}, U::AbstractGate{M}, N::Integer) where {K,M}

Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in quantum register `qreg` and affecting wires in tuple `itarget`.
"""
function controlled_circuit_gate(icntrl::NTuple{K, Union{Integer, Expr}}, qreg::QRegister, U::AbstractGate{M}, N::Integer; ccntrl::AbstractVector{<:Integer}=Int[]) where {K,M}
    !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    [controlled_circuit_gate(icntrl, i, U, N) for i in qreg.ind]
end

controlled_circuit_gate(icntrl::Union{<:Integer, Expr}, qreg::QRegister, U::AbstractGate{M}, N::Integer; ccntrl::AbstractVector{<:Integer}=Int[]) where {K,M} =
    controlled_circuit_gate((icntrl,), qreg, U, N)

"""
    controlled_circuit_gate(qreg1::QRegister, qreg2::QRegister, U::AbstractGate{M}, N::Integer) where {K,M}

Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in quantum register `qreg` and affecting wires in tuple `itarget`.
"""
function controlled_circuit_gate(reg::Register, qreg::QRegister, U::AbstractGate{M}, N::Integer; ccntrl::AbstractVector{<:Integer}=Int[]) where {K,M}
    !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    !any(reg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
    length(qreg.ind) == length(reg.ind) || error("Registers used must be of same length")
    [controlled_circuit_gate(i, j, U, N) for (i,j) in zip(reg.ind, qreg.ind)]
end


"""
    AbstractMoment{N}

Represents an intermediate state within a given circuit.
`N` is the overall number of quantum "wires" of the circuit.
"""
abstract type AbstractMoment{N} end

"""
    AbstractMoment{N}

Represents an intermediate state within a given circuit.
`N` is the overall number of quantum "wires" of the circuit.
"""
mutable struct Moment{N} <: AbstractMoment{N}
    gates::AbstractVector{<:AbstractCircuitGate{N}}

    @doc """
        Moment{N}(g::AbstractCircuitGate{N}) where {N}

    Create a `Moment{N}` object consisting of a single `CircuitGate{N}` object.
    """
    function Moment{N}(g::AbstractCircuitGate{N}) where {N}
        new([g])
    end

    @doc """
        Moment{N}(g::AbstractVector{<:AbstractCircuitGate{N}}) where {N}

    Create a `Moment{N}` object consisting of multiple `CircuitGate{N}` objects.
    """
    function Moment{N}(g::AbstractVector{<:AbstractCircuitGate{N}}) where {N}
        wires = Integer[]
        for gate in g
            length(intersect(wires, gate.iwire)) == 0 || error("Only gates on different wires are allowed in a Moment")
            append!(wires, gate.iwire)
        end
        new(g)
    end
end

"""
    Base.adjoint(m::Moment{N}) where {N}

Construct a `Moment{N}` object that is the adjoint of `m`.
"""
function Base.adjoint(m::Moment{N}) where {N}
    return Moment{N}(Base.adjoint.(reverse(m.gates)))
end

# make Moment iterable and indexable
function Base.getindex(m::Moment{N}, i::Integer) where {N}
    1 <= i <= length(m.gates) || throw(BoundsError(m, i))
    return m.gates[i]
end

function Base.iterate(m::Moment{N}, state=1) where {N}
    return state > length(m.gates) ? nothing : (m[state], state+1)
end

# implement methods required for iteration
function Base.firstindex(m::Moment{N}) where {N}
    return 1
end

function Base.lastindex(m::Moment{N}) where {N}
    return length(m.gates)
end

function Base.length(m::Moment{N}) where {N}
    return length(m.gates)
end

function Base.isapprox(m1::Moment{N}, m2::Moment{N}) where {N}
    length(m1) == length(m2) || return false
    for (g1, g2) in zip(m1, m2)
        if !(g1 ≈ g2)
            return false
        end
    end
    return true
end

function Base.reverse(m::Moment{N}) where {N}
    m.gates = reverse(m.gates)
    return m
end


"""
    rdm(N, iwire, ψ, χ)

Compute the reduced density matrix ``tr_B[|ψ⟩⟨χ|]``, where the trace runs over
the subsystem complementary to the qubits specified by `iwire`.
"""
function rdm(N::Integer, iwire::NTuple{M,<:Integer}, ψ::AbstractVector, χ::AbstractVector) where {M}
    M ≥ 1 || error("Need at least one wire to act on.")
    M ≤ N || error("Number of gate wires cannot be larger than total number of wires.")
    length(unique(iwire)) == M || error("Wire indices must be unique.")
    minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")
    maximum(iwire) ≤ N || error("Wire index cannot be larger than total number of wires.")

    # convert to array
    iwire = collect(iwire)
    # complementary wires
    iwcompl = setdiff(1:N, iwire)
    @assert length(iwire) + length(iwcompl) == N

    # TODO: support general "qudits"
    d = 2

    ρ = zeros(eltype(ψ), d^M, d^M)

    # Note: following the ordering convention of `kron` here, i.e.,
    # last qubit corresponds to fastest varying index
    strides = [d^(j-1) for j in 1:N]
    wstrides = strides[iwire]
    cstrides = strides[iwcompl]

    # TODO: optimize memory access pattern
    for kw in (N-M > 0 ? reverse.(cartesian_tuples(d, N-M)) : [Int[]])
        koffset = dot(collect(kw), strides[iwcompl])
        for (i, iw) in enumerate(reverse.(cartesian_tuples(d, M)))
            for (j, jw) in enumerate(reverse.(cartesian_tuples(d, M)))
                rowind = koffset + dot(collect(iw), strides[iwire]) + 1
                colind = koffset + dot(collect(jw), strides[iwire]) + 1
                ρ[i, j] += ψ[rowind] * conj(χ[colind])
            end
        end
    end

    return ρ
end

"""
    CircuitGateChain{N}

Chain of quantum circuit gates in a circuit of `N` qubits.
"""
mutable struct CircuitGateChain{N}
    moments::AbstractVector{<:AbstractMoment{N}}
    "classical registers"
    creg::AbstractVector{BitArray{1}}

    @doc """
        CircuitGateChain{N}(gates::AbstractVector{<:AbstractCircuitGate{N}}; creg::AbstractVector{BitArray}=BitArray[]) where {N}

    Chain of quantum circuit gates in a circuit of `N` qubits. Constructed from vector of `CircuitGate{N}` objects.
    """
    function CircuitGateChain{N}(gates::AbstractVector{<:AbstractCircuitGate{N}}) where {N}
        moments = map(Moment{N}, gates)
        new(moments, BitArray{1}[])
    end

    @doc """
        CircuitGateChain{N}(moments::AbstractVector{<:AbstractMoment{N}}; creg::AbstractVector{BitArray}=BitArray[]) where {N}

    Chain of quantum circuit gates in a circuit of `N` qubits, constructing from vector of `Moment{N}` objects.
    """
    function CircuitGateChain{N}(moments::AbstractVector{<:AbstractMoment{N}}; creg::Union{AbstractVector{BitArray},AbstractVector{Integer}}=BitArray[]) where {N}
        new(moments, BitArray{1}[])
    end

    @doc """
        CircuitGateChain(qregs::AbstractVector{QRegister}, cregs::AbstractVector{CRegister})

    Chain of quantum circuit gates in a circuit of `N` qubits, constructing from vectors quantum registers and classical registers.
    """
    function CircuitGateChain(qregs::AbstractVector{QRegister}, cregs::AbstractVector{CRegister})
        N = 1
        ccount = 1
        creg = BitArray{1}[]
        for reg in qregs
            reg.ind .= [N:N+reg.n-1...]
            N += reg.n
        end
        for reg in cregs
            reg.ind .= [-ccount:-1:-ccount-length(reg.n)+1...]
            push!(creg, reg.n)
            ccount += length(reg.n)
        end
        cgc = new{N-1}(Moment{N-1}[], creg)
    end

    @doc """
        CircuitGateChain(qregs::AbstractVector{QRegister})

    Chain of quantum circuit gates in a circuit of `N` qubits, constructing from vectors quantum registers.
    """
    function CircuitGateChain(qregs::AbstractVector{QRegister})
        N = 1
        creg = BitArray{1}[]
        for reg in qregs
            reg.ind .= [N:N+reg.n-1...]
            N += reg.n
        end
        cgc = new{N-1}(Moment{N-1}[], creg)
    end

    @doc """
        CircuitGateChain{N}(cregs::AbstractVector{CRegister})

    Chain of quantum circuit gates in a circuit of `N` qubits, constructing with vector of classical registers.
    """
    function CircuitGateChain{N}(cregs::AbstractVector{CRegister}) where {N}
        ccount = 1
        creg = BitArray{1}[]
        for reg in cregs
            reg.ind .= [-ccount:-1:-ccount-length(reg.n)+1...]
            push!(creg, reg.n)
            ccount += length(reg.n)
        end
        cgc = new{N}(Moment{N}[], creg)
    end

    @doc """
        CircuitGateChain{N}()

    Chain of quantum circuit gates in a circuit of `N` qubits.
    """
    function CircuitGateChain{N}() where {N}
        new{N}(Moment{N}[], BitArray{1}[])
    end

end

"""
    add_creg!(cgc::CircuitGateChain{N}, reg::CRegister) where {N}

Add classical register to CircuitGateChain `cgc`
"""
function add_creg!(cgc::CircuitGateChain{N}, reg::CRegister) where {N}
    all(reg.ind .== 0) || error("CRegister is already used in another CircuitGateChain object")
    if !isempty(cgc.creg)
        ccount = length(collect(Iterators.flatten(reverse.(cgc.creg)))) + 1
    else
        ccount = 1
    end
    reg.ind .= [-ccount:-1:-ccount-length(reg.n)+1...]

    push!(cgc.creg, reg.n)
end


# unitary matrix representation of a sequence of circuit gates
"""
    matrix(cgc::CircuitGateChain{N}) where {N}

returns matrix representation of a `CircuitGateChain{N}` object that can be applied to a state vector of `N` qubits.
"""
function matrix(cgc::CircuitGateChain{N}) where {N}
    gates = AbstractCircuitGate{N}[]
    for moment in cgc
        append!(gates, moment.gates)
    end
    prod(Tuple(matrix(g) for g in reverse(gates)))
end

"""
    Base.adjoint(cgc::CircuitGateChain{N}) where {N}

Construct a `CircuitGateChain{N}` object that is the adjoint of `cgc`.
"""
function Base.adjoint(cgc::CircuitGateChain{N}) where {N}
    return CircuitGateChain{N}(Base.adjoint.(reverse(cgc.moments)))
end

# make CircuitGateChain iterable and indexable
function Base.getindex(cgc::CircuitGateChain{N}, i::Integer) where {N}
    1 <= i <= length(cgc.moments) || throw(BoundsError(cgc, i))
    return cgc.moments[i]
end

function Base.setindex!(cgc::CircuitGateChain{N}, m::Moment{N}, i::Integer) where {N}
    1 <= i <= length(cgc.moments) || throw(BoundsError(cgc, i))
    cgc.moments[i] = m
end

function Base.iterate(cgc::CircuitGateChain{N}, state=1) where {N}
    return state > length(cgc.moments) ? nothing : (cgc[state], state+1)
end

# implement methods required for iteration
function Base.firstindex(cgc::CircuitGateChain{N}) where {N}
    return 1
end

function Base.lastindex(cgc::CircuitGateChain{N}) where {N}
    return length(cgc.moments)
end

function Base.length(cgc::CircuitGateChain{N}) where {N}
    return length(cgc.moments)
end

Base.size(cgc::CircuitGateChain{N}) where {N} = N

function Base.:*(cgc1::CircuitGateChain{N}, cgc2::CircuitGateChain{N}) where {N}
    append!(cgc1.moments, cgc2.moments)
    creg_length = length(cgc1.creg) - length(cgc2.creg)
    if creg_length > 0
        append!(cgc2.creg, fill(0, (creg_length,)))
    elseif creg_length < 0
        append!(cgc1.creg, fill(0, (abs(creg_length)),))
    end
    cgc1.creg = cgc1.creg .| cgc2.creg
    return cgc1
end

function Base.reverse(cgc::CircuitGateChain{N}) where {N}
    for i in length(cgc)
        cgc[i] = reverse(cgc[i])
    end
    cgc.moments = reverse(cgc.moments)
    return cgc
end


function (cgc::CircuitGateChain{N})(g::CircuitGate{M,N,G}) where {M,N,G<:AbstractGate}
    max_creg = sum(length.(cgc.creg))
    if g.ccntrl isa Vector{Integer}
        for reg in g.ccntrl
            reg <= max_creg || error("Attempt to access classical bit " * string(reg) * " in CircuitGateChain with " * string(max_creg) * " classical bits")
        end
    end

    push!(cgc.moments, Moment{N}(g))
end

function (cgc::CircuitGateChain{N})(g::CircuitGate{M,N,G}, max_creg::Int) where {M,N,G<:AbstractGate}
    if g.ccntrl isa Vector{Integer}
        for reg in g.ccntrl
            reg <= max_creg || error("Attempt to access classical bit " * string(reg) * " in CircuitGateChain with " * string(max_creg) * " classical bits")
        end
    end

    push!(cgc.moments, Moment{N}(g))
end

function (cgc::CircuitGateChain{N})(g::Array{<:CircuitGate{M,N,G} where M where G,1}) where {N}
    max_creg = sum(length.(cgc.creg))
    for cg in g
        cgc(cg, max_creg)
    end
end

function (cgc::CircuitGateChain{N})(g::Array{<:Any,1}) where {N}
    max_creg = sum(length.(cgc.creg))
    cgc(g, max_creg)
end

function (cgc::CircuitGateChain{N})(g::Array{<:Any,1}, max_creg::Int) where {N}
    for cg in g
        if cg isa AbstractVector
            cgc.(cg, max_creg)
            continue
        end
        cgc(cg, max_creg)
    end
end

"""
    MeasurementOps{N}

Pairwise commuting measurement operators (Hermitian matrices) for circuit of size `N`.
"""
struct MeasurementOps{N}
    mops::AbstractVector{<:AbstractMatrix}
    cgs::AbstractVector{<:CircuitGate}

    @doc """
        MeasurementOps{N}(mop::AbstractMatrix) where {N}

    constructs a `MeasurementOps{N}` object for circuit of size `N` from single ``2^{N} \\times 2^{N}`` matrix.
    """
    function MeasurementOps{N}(mop::AbstractMatrix) where {N}
        # TODO: support general "qudits"
        d = 2
        # consistency checks
        size(mop) == (d^N, d^N) || error("Measurement operator must be a 2^N × 2^N matrix.")
        mop ≈ Base.adjoint(mop) || error("Measurement operator must be Hermitian.")
        new([mop])
    end

    @doc """
        MeasurementOps{N}(mops::AbstractVector{<:AbstractMatrix}) where {N}

    constructs a `MeasurementOps{N}` object for a circuit of size `N` from a vector of ``2^{N} \\times 2^{N}`` matrices.
    """
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

    @doc """
        MeasurementOps{N}(cg::CircuitGate) where {N}

    constructs a `MeasurementOps{N}` object for a circuit of size `N` from a single `CircuitGate{N}` object.
    """
    function MeasurementOps{N}(cg::CircuitGate) where {N}
        # TODO: support general "qudits"
        d = 2
        # consistency checks
        mop = matrix(cg)
        size(mop) == (d^N, d^N) || error("Measurement operator must be a 2^N × 2^N matrix.")
        mop ≈ Base.adjoint(mop) || error("Measurement operator must be Hermitian.")
        new([mop], [cg])
    end

    @doc """
        MeasurementOps{N}(cgs::AbstractVector{<:CircuitGate}) where {N}

    constructs a `MeasurementOps{N}` object for a circuit of size `N` from a vector of `CircuitGate{N}` objects.
    """
    function MeasurementOps{N}(cgs::AbstractVector{<:CircuitGate}) where {N}
        # TODO: support general "qudits"
        d = 2
        # consistency checks
        mops = AbstractMatrix[]
        for cg in cgs
            push!(mops, matrix(cg))
        end
        for m in mops
            size(m) == (d^N, d^N) || error("Measurement operator must be a 2^N × 2^N matrix.")
            m ≈ Base.adjoint(m) || error("Measurement operator must be Hermitian.")
            for n in mops
                norm(comm(m, n))/d^N < 1e-13 || error("Measurement operators must pairwise commute.")
            end
        end
        new(mops, cgs)
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

# make CircuitGateChain iterable and indexable
function Base.getindex(c::Circuit{N}, i::Integer) where {N}
    1 <= i <= length(c.cgc.moments) || throw(BoundsError(c, i))
    return c.cgc.moments[i]
end

function Base.iterate(c::Circuit{N}, state=1) where {N}
    return state > length(c.cgc.moments) ? nothing : (c.cgc[state], state+1)
end

# implement methods required for iteration
function Base.firstindex(c::Circuit{N}) where {N}
    return 1
end

function Base.lastindex(c::Circuit{N}) where {N}
    return length(c.cgc.moments)
end

function distribution(c::Circuit{N}, ψ::AbstractVector) where {N}
    return apply(c.cgc, ψ)
end
