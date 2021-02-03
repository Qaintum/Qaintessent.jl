using LinearAlgebra

"""
    cartesian_tuples(d, N)

Enumerate cartesian tuples.

# Example
```julia-repl
julia> cartesian_tuples(2, 3)
2×2×2 Array{Tuple{Int64,Int64,Int64},3}:
[:, :, 1] =
 (0, 0, 0)  (0, 1, 0)
 (1, 0, 0)  (1, 1, 0)

[:, :, 2] =
 (0, 0, 1)  (0, 1, 1)
 (1, 0, 1)  (1, 1, 1)
```
"""
cartesian_tuples(d::Integer, N::Integer) =
    cartesian_tuples(d, Val(N))

cartesian_tuples(d::Integer, ::Val{N}) where {N} =
    Tuple.(CartesianIndices(NTuple{N, UnitRange{Int64}}(fill(0:d - 1, N))))


"""
    comm(A, B)

Matrix commutator [A, B].
"""
comm(A::AbstractMatrix, B::AbstractMatrix) = A * B - B * A


"""
    pauli_vector(x, y, z)

Assemble the "Pauli vector" matrix.
"""
pauli_vector(x, y, z) = ComplexF64[z x - im * y; x + im * y -z]


"""
    binary_rep(x::Integer, M::Integer)

Return binary representation of Integer `x` for `M` bits with least significant bit first.
"""
function binary_rep(x::Integer, M::Integer)
    m = BitArray(undef, M)
    for i in 1:M
        m[i] = x & 1
        x = x >> 1
    end
    m
end

"""
    binary_rep!(x::Integer, M::Integer)

Return binary representation of Integer `x` for `M` bits with least significant bit first.
"""
function binary_rep!(m::BitArray{1}, x::Integer, M::Integer)
    for i in 1:M
        m[i] = x & 1
        x = x >> 1
    end
    m
end


"""
    intlog2(x::Integer)

returns log 2 with integer inputs
"""
function intlog2(x::Integer)
    x == 0 && error("Logarithm of 0 is undefined")
    x < 0 && error("Logarithm of negative number is undefined")
    ret::Int = 0
    while x > 1
        x = x >> 1
        ret += 1
    end
    ret
end


"""
    sliced_index(idx::Tuple, targetwires::Tuple, N::Int)

Construct sliced index via target map (utility function)
"""
function sliced_index(idx::Tuple, targetwires::Tuple, N::Int)
    islice = Vector{Any}(fill(Colon(), N))
    M = length(targetwires)
    for k in 1:M
        islice[targetwires[k]] = idx[k] + 1
    end
    return islice
end


"""
    gramm_schmidt(a::Array{ComplexF64})

returns orthogonalized complex vectors
"""
function gramm_schmidt!(a::Array{ComplexF64})
    num_vectors = size(a,2)
    a[:,end] = a[:,end] ./ norm(a[:,end])
    for i in num_vectors-1:-1:1
        while !(norm(dot(a[:,i+1],a[:,i])) < 1e-16)
            a[:, i] -= conj(dot(a[:, i],a[:, i+1]) / dot(a[:, i], a[:, i])) * a[:, i+1]
            a[:, i] = a[:, i] ./ norm(a[:, i])
        end
    end
    a
end

"""
    gramm_schmidt(a::Array{Float64})

returns orthogonalized real vectors
"""
function gramm_schmidt!(a::Array{Float64})
    num_vectors = size(a,2)
    a[:,end] = a[:,end] ./ norm(a[:,end])
    for i in num_vectors-1:-1:1
        while !(norm(dot(a[:,i+1],a[:,i])) < 1e-16)
            a[:, i] -= dot(a[:, i],a[:, i+1] / dot(a[:, i], a[:, i])) * a[:, i+1]
            a[:, i] = a[:, i] ./ norm(a[:, i])
        end
    end
    a
end