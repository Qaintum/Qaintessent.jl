using LinearAlgebra
using Memoize

"""
    pauli_vector(x, y, z)

Assemble the "Pauli vector" matrix.
"""
pauli_vector(x, y, z) = ComplexQ[z x - im * y; x + im * y -z]


"""
    binary_digits(x::Integer, M::Integer)

Return `M` bits of the binary representation of the integer `x` as array (least significant bit at first index).
"""
function binary_digits(M::Integer, x::Integer)
    binary_digits!(BitArray(undef, M), x)
end

"""
    binary_digits!(m::BitArray{1}, x::Integer)

Fill `m` with the binary representation of the integer `x` (least significant bit at first index).
"""
function binary_digits!(m::BitArray{1}, x::Integer=0)
    for i in 1:length(m)
        m[i] = x & 1
        x >>= 1
    end
    return m
end

"""
    binary_to_int(m::BitArray{1}, s=0)

converts BitArray `m` into Integer
"""
function binary_to_int(m::BitArray{1}, x=0)
    v = 1
    for i in Base.view(m, :)
        x += v*i
        v <<= 1
    end 
    x
end


"""
    quaternary_digits(x::Integer, M::Integer)

Return `M` base-4 digits of the quaternary representation of the integer `x` as array (least significant digit at first index).
"""
function quaternary_digits(M::Integer, x::Integer)
    quaternary_digits!(Vector{Int}(undef, M), x)
end

"""
    quaternary_digits!(m::Vector{Int}, x::Integer)

Fill `m` with the quaternary representation of the integer `x` (least significant digit at first index).
"""
function quaternary_digits!(m::Vector{Int}, x::Integer)
    for i in 1:length(m)
        m[i] = x & 3
        x >>= 2
    end
    return m
end


"""
    intlog2(x::Integer)

Compute integer base-2 logarithm. Rounds to floor if x is not a power of 2.
"""
function intlog2(x::Integer)
    x == 0 && error("Logarithm of 0 is undefined")
    x < 0 && error("Logarithm of negative number is undefined")
    ret::Int = 0
    while x > 1
        x >>= 1
        ret += 1
    end
    ret
end


"""
    sliced_index(idx, targetwires::Tuple, N::Int)

Construct sliced index via target map (utility function).
"""
function sliced_index(idx, targetwires::Tuple, N::Int)
    islice = Vector{Any}(fill(Colon(), N))
    for k in 1:length(targetwires)
        islice[targetwires[k]] = idx[k] + 1
    end
    return islice
end


"""
    gramm_schmidt(a::Array{<:Real})

Compute orthogonalized real vectors.
"""
function gramm_schmidt!(a::Array{<:Real})
    num_vectors = size(a,2)
    
    a[:, end] = a[:, end] ./ norm(a[:, end])

    for i in num_vectors-1:-1:1        
        for j in i+1:num_vectors
            a[:, i] = a[:, i] - (dot(a[:, i], a[:, j]) / dot(a[:, j], a[:, j])) * a[:, j]
        end
        a[:, i] = a[:, i] ./ norm(a[:, i])
    end

    return a
end

"""
    gramm_schmidt(a::Array{<:Complex})

Compute orthogonalized complex vectors.
"""
function gramm_schmidt!(a::Array{<:Complex})
    num_vectors = size(a,2)

    a[:, end] = a[:, end] ./ norm(a[:, end])

    for i in num_vectors-1:-1:1        
        for j in num_vectors:-1:i+1
            a[:, i] = a[:, i] - conj(dot(a[:, i], a[:, j]) / dot(a[:, j], a[:, j])) * a[:, j]
        end
        a[:, i] = a[:, i] ./ norm(a[:, i])
    end
    return a
end

"""
    rdm(N, iwire, ψ, χ, d=2)

Compute the reduced density matrix ``tr_B[|ψ⟩⟨χ|]``, where the trace runs over
the subsystem complementary to the qubits specified by `iwire`.
"""
@views function rdm(N::Integer, iwire::NTuple{M,<:Integer}, ψ::Vector{T}, χ::Vector{T}, d::Int=2) where {M,T}
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
    # # Note: following the ordering convention of `kron` here, i.e.,
    # # last qubit corresponds to fastest varying index
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
                    @inbounds ρ[i, j] += ψ[rowind] * conj(χ[colind])
                end
            end
        end
    elseif d == 4
        iw = Vector{Int}(undef, M)
        jw = Vector{Int}(undef, M)
        kw = Vector{Int}(undef, N - M)
        # TODO: optimize memory access pattern
        for k in 0:4^(N-M)-1
            quaternary_digits!(kw, k)
            koffset = dot(kw, cstrides)
            for i in 1:4^M
                quaternary_digits!(iw, i - 1)
                for j in 1:4^M
                    quaternary_digits!(jw, j - 1)
                    rowind = koffset + dot(iw, wstrides) + 1
                    colind = koffset + dot(jw, wstrides) + 1
                    @inbounds ρ[i, j] += ψ[rowind] * conj(χ[colind])
                end
            end
        end
    else
        error("d = $d not supported yet.")
    end

    return ρ
end

data(::AbstractMatrix) = nothing