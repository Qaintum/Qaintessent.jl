using LinearAlgebra


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
    binary_digits(x::Integer, M::Integer)

Return `M` bits of the binary representation of the integer `x` as array (least significant bit at first index).
"""
function binary_digits(x::Integer, M::Integer)
    binary_digits!(BitArray(undef, M), x)
end

"""
    binary_digits!(m::BitArray{1}, x::Integer)

Fill `m` with the binary representation of the integer `x` (least significant bit at first index).
"""
function binary_digits!(m::BitArray{1}, x::Integer)
    for i in 1:length(m)
        m[i] = x & 1
        x >>= 1
    end
    return m
end


"""
    quaternary_digits(x::Integer, M::Integer)

Return `M` base-4 digits of the quaternary representation of the integer `x` as array (least significant digit at first index).
"""
function quaternary_digits(x::Integer, M::Integer)
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

Compute integar base-2 logarithm.
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
    gramm_schmidt(a::Array{Float64})

Compute orthogonalized real vectors.
"""
function gramm_schmidt!(a::Array{Float64})
    num_vectors = size(a,2)
    a[:, end] = a[:, end] ./ norm(a[:, end])
    for i in num_vectors-1:-1:1
        while norm(dot(a[:, i+1], a[:, i])) > 1e-16
            a[:, i] -= dot(a[:, i], a[:, i+1] / dot(a[:, i], a[:, i])) * a[:, i+1]
            a[:, i] = a[:, i] ./ norm(a[:, i])
        end
    end
    return a
end

"""
    gramm_schmidt(a::Array{ComplexF64})

Compute orthogonalized complex vectors.
"""
function gramm_schmidt!(a::Array{ComplexF64})
    num_vectors = size(a, 2)
    a[:, end] = a[:, end] ./ norm(a[:, end])
    for i in num_vectors-1:-1:1
        while norm(dot(a[:, i+1], a[:, i])) > 1e-16
            a[:, i] -= conj(dot(a[:, i], a[:, i+1]) / dot(a[:, i], a[:, i])) * a[:, i+1]
            a[:, i] = a[:, i] ./ norm(a[:, i])
        end
    end
    return a
end
