using Memoize
"""
    Statevector

Represents a pure state.
"""
struct Statevector <: AbstractArray{ComplexF64,1}
    state::Vector{ComplexF64}
    vec::Vector{ComplexF64}
    perm::Vector{Int}
    N::Int

    @doc """
        Statevector(v::Vector{ComplexF64})

    Create a `Statevector` object representing a pure quantum state.
    """
    function Statevector(v::Vector{ComplexF64})
        ispow2(length(v)) || error("Statevector must have length which is a power of 2")
        N = intlog2(length(v))
        new(v, zeros(ComplexF64, 2^N), collect(1:2^N), N)
    end
    @doc """
    Statevector(v::Vector{ComplexF64})

    Create a `Statevector` object representing a pure quantum state with 'N' qubits.
    """
    function Statevector(N::Int)
        new(zeros(ComplexF64, 2^N), zeros(ComplexF64, 2^N), collect(1:2^N), N)
    end
end



Base.setindex!(ψ::Statevector, v, i::Int64) = setindex!(ψ.state, v, i)
Base.getindex(ψ::Statevector, i::Int) = getindex(ψ.state, i)
Base.size(ψ::Statevector) = size(ψ.state)
Base.length(ψ::Statevector) = length(ψ.state)
Base.IndexStyle(ψ::Type{Statevector}) = IndexLinear()

"""
    rdm(N, iwire, ψ, χ, d=2)

Compute the reduced density matrix ``tr_B[|ψ⟩⟨χ|]``, where the trace runs over
the subsystem complementary to the qubits specified by `iwire`.
"""
function rdm(N::Integer, iwire::NTuple{M,<:Integer}, ψ::Statevector, χ::Statevector, d::Int=2) where {M}
    M ≥ 1 || error("Need at least one wire to act on.")
    M ≤ N || error("Number of gate wires cannot be larger than total number of wires.")
    length(unique(iwire)) == M || error("Wire indices must be unique.")
    minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")
    maximum(iwire) ≤ N || error("Wire index cannot be larger than total number of wires.")
    length(ψ) == d^N || error("Input vector 'ψ' has wrong length.")
    length(χ) == d^N || error("Input vector 'χ' has wrong length.")

    # convert to array
    iwire = collect(iwire)
    # complementary wires
    iwcompl = setdiff(1:N, iwire)
    @assert length(iwire) + length(iwcompl) == N

    ρ = zeros(eltype(ψ), d^M, d^M)

    # Note: following the ordering convention of `kron` here, i.e.,
    # last qubit corresponds to fastest varying index
    strides = [d^(j-1) for j in 1:N]
    wstrides = strides[iwire]
    cstrides = strides[iwcompl]

    if d == 2
        iw = BitArray{1}(undef, M)
        jw = BitArray{1}(undef, M)
        kw = BitArray{1}(undef, N - M)
        # TODO: optimize memory access pattern
        for k in 0:2^(N-M)-1
            binary_digits!(kw, k)
            koffset = dot(kw, cstrides)
            for i in 1:2^M
                binary_digits!(iw, i - 1)
                for j in 1:2^M
                    binary_digits!(jw, j - 1)
                    rowind = koffset + dot(iw, wstrides) + 1
                    colind = koffset + dot(jw, wstrides) + 1
                    ρ[i, j] += ψ[rowind] * conj(χ[colind])
                end
            end
        end
    else
        error("d = $d not supported yet.")
    end

    return ρ
end