"""
AbstractCircuitGate

Abstract unitary quantum circuit gate. `N` is the overall number of quantum "wires" of the circuit.
"""
abstract type AbstractCircuitGate end

"""
CircuitGate{M,G} <: AbstractCircuitGate

Unitary quantum circuit gate. `M` is the number of wires affected by the CircuitGate, `N` is the overall number of quantum "wires" of the circuit, `G` is the basic gate used to construct the CircuitGate.
"""
struct CircuitGate{M,G} <: AbstractCircuitGate
    "ordered wire indices which this gate acts on"
    iwire::NTuple{M,Int64}
    "abstract gate"
    gate::G

    @doc """
        CircuitGate{M,G}(iwire::NTuple{M,<:Integer}, gate::G) where {M,G}

    Creates a `CircuitGate{M,G}` object. `M` is the number of wires affected by the CircuitGate, `N` is the overall number of quantum "wires" of the circuit, `G` is the basic gate used to construct the CircuitGate.
    """
    function CircuitGate{M,G}(iwire::NTuple{M,<:Integer}, gate::G) where {M,G <: AbstractGate}
        M == wires(gate) || error("$G affects $(wires(gate)) wires but $M wires, $iwire, were passed.")
        M ≥ 1 || error("Need at least one wire to act on.")
        length(unique(iwire)) == M || error("Wire indices must be unique.")
        minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")

        new{M,G}(iwire, gate)
    end
end

"""
CircuitGate(iwire::NTuple{M,<:Integer}, gate::AbstractGate, N) where {M}

Creates a `CircuitGate{M,N,G}` object. `M` is the number of wires affected by the CircuitGate, `N` is the overall number of quantum "wires" of the circuit, `G` is the basic gate used to construct the CircuitGate.
"""
function CircuitGate(iwire::NTuple{M,<:Integer}, gate::G) where {M,G <:AbstractGate}
    CircuitGate{M,G}(iwire, gate)
end

data(cg::CircuitGate) = data(cg.gate)

"""
Base.isapprox(cg1::CircuitGate{M,G}, cg2::CircuitGate{M,G})

Compares two circuit gates of basic type `G`.
"""
function Base.isapprox(cg1::CircuitGate{M,G}, cg2::CircuitGate{M,G}) where {M,G}
# for parametric gates, return true only if parameters are approximately equal
    if cg1.iwire != cg2.iwire
        return false
    end
    return cg1.gate ≈ cg2.gate
end

Base.isapprox(cg1::CircuitGate{M,G}, cg2::CircuitGate{M,H}) where {M,G,H} = false

LinearAlgebra.ishermitian(cg::CircuitGate) = LinearAlgebra.ishermitian(cg.gate)


"""
matrix(cg::CircuitGate{M,G}) where {M,G<:AbstractGate}

returns matrix representation of a `CircuitGate{M,G}` object that can applied to a state vector of `N` qubits.
"""
function matrix(cg::CircuitGate{M,G}, N::Integer=0) where {M,G <: AbstractGate}
# convert to array
    iwire = Int64[i for i in cg.iwire]
    if N == 0
        N = maximum([iwire' M])
    else
        N >= maximum(iwire) || error("CircuitGate applied to iwires, $iwire. Input circuit size `N` is $N")
    end
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
    strides = Int64[d^(j - 1) for j in 1:N]
    wstrides = strides[iwire]
    cstrides = strides[iwcompl]
    b = BitArray(undef, M)

    values = Vector{ComplexF64}(undef, length(gmat.nzval)*d^(N-M))
    value_size = length(gmat.nzval)
    rowind = Vector{Int64}(undef, length(values))
    colind = Vector{Int64}(undef, length(values))

    for j in 1:length(gmat.colptr)-1
        index = gmat.colptr[j]:gmat.colptr[j+1]-1
        colvalue = dot(binary_rep!(b, j-1, M), wstrides) + 1
        for i in index
            @inbounds rowind[i] = dot(binary_rep!(b, gmat.rowval[i] - 1, M), wstrides) + 1
            @inbounds colind[i] = colvalue
        end
    end

    r = @view rowind[1:value_size]
    c = @view colind[1:value_size]
    @inbounds values[1:value_size] = gmat.nzval
    v = @view values[1:value_size]

    b = BitArray(undef, N-M)
    for kw in 2:d^(N-M)
        koffset = dot(binary_rep!(b, kw-1, N-M), cstrides)
        @inbounds rowind[value_size*(kw-1)+1:value_size*(kw)] = r .+ koffset
        @inbounds colind[value_size*(kw-1)+1:value_size*(kw)] = c .+ koffset
        @inbounds values[value_size*(kw-1)+1:value_size*(kw)] = v
    end

    return sparse(rowind, colind, values, d^N, d^N)
end

"""
Base.adjoint(cg::CircuitGate{M,G})

Construct a `CircuitGate{M,H}` object where `H` is the adjoint of `AbstractGate` `G`
"""
function Base.adjoint(cg::CircuitGate{M,G}) where {M,G}
    adj_gate = Base.adjoint(cg.gate)
    CircuitGate{M,typeof(adj_gate)}(cg.iwire, adj_gate)
end


"""
single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate{1}, N::Integer)

Construct a `CircuitGate{1,G}` object of basic gate type `gate` affecting wire `iwire`.
"""
single_qubit_circuit_gate(iwire::Integer, gate::AbstractGate) =
CircuitGate((iwire,), gate)

# """
#     single_qubit_circuit_gate(qreg::QRegister, gate::AbstractGate{1}, N::Integer)

# Construct a `CircuitGate{1,N,G}` objects of basic gate type `gate` affecting QRegister `qreg`
# """
# function single_qubit_circuit_gate(qreg::QRegister, gate::AbstractGate{1}, N::Integer)
#     !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     single_qubit_circuit_gate.(qreg.ind, (gate,), (N,))
# end

function cg(iwire::Integer, gate::AbstractGate, control::NTuple{M,Integer}=NTuple{0,Integer}()) where {M}
    if isempty(control)
        return CircuitGate((iwire,), gate)
    else
        return CircuitGate((iwire, control...), ControlledGate(gate, M))
    end
end

function cg(iwire::Integer, gate::AbstractGate, control::Integer...)
    cg(iwire, gate, control)
end

"""
two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate{2}, N::Integer)

Construct a `CircuitGate{2,G}` object of basic gate type `gate` affecting wires `iwire1` and `iwire2`.
"""
two_qubit_circuit_gate(iwire1::Integer, iwire2::Integer, gate::AbstractGate) =
CircuitGate((iwire1, iwire2), gate)

function cg(iwire1::Integer, iwire2::Integer, gate::AbstractGate, control::NTuple{M,Integer}=NTuple{0,Integer}()) where {M}
    if isempty(control)
        return CircuitGate((iwire1, iwire2), gate)
    else
        C = length(control)
        return CircuitGate((iwire1, iwire2, control...), ControlledGate(gate, C))
    end
end

function cg(iwire1::Integer, iwire2::Integer, gate::AbstractGate, control::Integer...)
    cg(iwire1, iwire2, gate, control)
end

# """
#     two_qubit_circuit_gate(qreg::QRegister, iwire2::Integer, gate::AbstractGate{2}, N::Integer)

# Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting QRegister `qreg` and wire `iwire2`.
# """
# function two_qubit_circuit_gate(qreg::QRegister, iwire2::Integer, gate::AbstractGate, N::Integer)
#     !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     two_qubit_circuit_gate.(qreg.ind, (iwire2,), (gate,), (N,))
# end

# """
#     two_qubit_circuit_gate(iwire1::Integer, qreg::QRegister, gate::AbstractGate{2}, N::Integer)

# Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting QRegister `qreg` and wire `iwire1`.
# """
# function two_qubit_circuit_gate(iwire1::Integer, qreg::QRegister, gate::AbstractGate, N::Integer)
#     !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     two_qubit_circuit_gate.((iwire1,), qreg.ind, (gate,), (N,))
# end

# """
#     two_qubit_circuit_gate(qreg1::QRegister, qreg1::QRegister, gate::AbstractGate{2}, N::Integer)

# Construct a `CircuitGate{2,N,G}` object of basic gate type `gate` affecting QRegister `qreg1` and ``qreg2`.
# """
# function two_qubit_circuit_gate(qreg1::QRegister, qreg2::QRegister, gate::AbstractGate{2}, N::Integer)
#     !any(qreg1.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     !any(qreg2.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     length(qreg1.ind) == length(qreg2.ind) || error("Only able to apply CircuitGate to two QRegister objects if registers are of same length")
#     two_qubit_circuit_gate.(qreg1.ind, qreg2.ind, (gate,), (N,))
# end


# single control and target wire
"""
controlled_circuit_gate(itarget::Integer, icntrl::Union{Integer, Expr}, U::AbstractGate{1}, N::Integer)

Construct a `CircuitGate{2,N,G}` object of basic gate type `U` controlled by wire or Expr `icntrl` and affecting wire `itarget`.
"""
controlled_circuit_gate(itarget::Integer, icntrl::Integer, U::AbstractGate) =
controlled_circuit_gate((itarget,), (icntrl,), U)

# single control wire
"""
controlled_circuit_gate(itarget::NTuple{M,<:Integer}, icntrl::Integer, U::AbstractGate{M}, N::Integer) where {M}

Construct a `CircuitGate{M+1,N,G}` object of basic gate type `U` controlled by wire or Expr `icntrl` and affecting wires in tuple `itarget`.
"""
controlled_circuit_gate(itarget::NTuple{M,<:Integer}, icntrl::Integer, U::AbstractGate) where {M} =
controlled_circuit_gate(itarget, (icntrl,), U)

# single target wire
"""
controlled_circuit_gate(itarget::Integer, icntrl::NTuple{K, Union{Int, Expr}}, U::AbstractGate{1}, N::Integer)  where {K}

Construct a `CircuitGate{K+1,N,G}` object of basic gate type `U` controlled by wires or Expr in tuple `icntrl` and affecting wire `itarget`.
"""
controlled_circuit_gate(itarget::Integer, icntrl::NTuple{K,Integer}, U::AbstractGate)  where {K} =
controlled_circuit_gate((itarget,), icntrl, U)

"""
controlled_circuit_gate(itarget::NTuple{M, <:Integer}, icntrl::NTuple{K, <:Union{Integer, Expr}}, U::AbstractGate{M}, N::Integer) where {K,M}

Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in tuple `icntrl` and affecting wires in tuple `itarget`.
"""
function controlled_circuit_gate(itarget::NTuple{M,<:Integer}, icntrl::NTuple{K,Integer}, U::AbstractGate) where {K,M}
    # # consistency checks
    C = length(icntrl)
    # C + M ≤ N || error("Number of control and target wires must be smaller than overall number of wires.")
    # length(intersect(itarget, icntrl)) == 0 || error("Control and target wires must be disjoint.")
    CircuitGate((itarget..., icntrl...), ControlledGate(U, C), N)
end

wires(cg::CircuitGate) = cg.iwire

# # control with QRegisters
# """
#     controlled_circuit_gate(itarget::NTuple{M, <:Integer}, reg::Register, U::AbstractGate{M}, N::Integer, ccntrl::AbstractVector{Integer}=Int[]) where {M}

# Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in register `reg` and affecting wires in tuple `itarget`.
# """
# function controlled_circuit_gate(itarget::NTuple{M,<:Integer}, reg::Register, U::AbstractGate{M}, N::Integer, ccntrl::AbstractVector{<:Integer}=Int[]) where {M}
#     !any(reg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     controlled_circuit_gate.((itarget,), reg.ind, (U,), (N,))
# end

# controlled_circuit_gate(itarget::Integer, reg::Register, U::AbstractGate{M}, N::Integer, ccntrl::AbstractVector{<:Integer}=Int[]) where {M} =
#     controlled_circuit_gate((itarget,), reg, U, N)

# """
#     controlled_circuit_gate(qreg::QRegister, icntrl::NTuple{M,<:Integer}, U::AbstractGate{M}, N::Integer) where {M}

# Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` affecting wires in quantum register `qreg` and controlled by wires in tuple `icntrl`.
# """
# function controlled_circuit_gate(qreg::QRegister, icntrl::NTuple{K,Union{Integer,Expr}}, U::AbstractGate{M}, N::Integer, ccntrl::AbstractVector{<:Integer}=Int[]) where {K,M}
#     !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     controlled_circuit_gate.(qreg.ind, (icntrl,), (U,), (N,))
# end

# controlled_circuit_gate(qreg::QRegister, icntrl::Union{<:Integer,Expr}, U::AbstractGate{M}, N::Integer, ccntrl::AbstractVector{<:Integer}=Int[]) where {M} =
#     controlled_circuit_gate(qreg, (icntrl,), U, N)

# """
#     controlled_circuit_gate(qreg1::QRegister, qreg2::QRegister, U::AbstractGate{M}, N::Integer) where {K,M}

# Construct a `CircuitGate{M+K,N,G}` object of basic gate type `U` controlled by wires in quantum register `qreg` and affecting wires in tuple `itarget`.
# """
# function controlled_circuit_gate(reg::Register, qreg::QRegister, U::AbstractGate{M}, N::Integer, ccntrl::AbstractVector{<:Integer}=Int[]) where {M}
#     !any(qreg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     !any(reg.ind .== 0) || error("Register object has yet to be used in a CircuitGateChain object")
#     length(qreg.ind) == length(reg.ind) || error("Registers used must be of same length")
#     controlled_circuit_gate.(reg.ind, qreg.ind, (U,), (N,))
# end