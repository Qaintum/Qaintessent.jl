
"""
    DensityMatrix

Density matrix, represented with respect to identity and Pauli basis (σ_j/2) by storing
4 real coefficients for each qubit, i.e., a real vector of length `4^N` for `N` qubits.
"""
mutable struct DensityMatrix
    "coefficients with respect to identity and Pauli basis"
    v::Vector{FloatQ}
    "number of qubits"
    N::Int
    "scratch space for memory efficient permutation"
    scratch::Vector{FloatQ}

    @doc """
        DensityMatrix(v::AbstractVector{<:Real}, N::Integer)

    Density matrix representation given Pauli-vector dotted with with vector of coefficients
    """
    function DensityMatrix(v::AbstractVector{<:Real}, N::Integer)
        length(v) == 4^N || error("Expected length of coefficient vector for density matrix is `4^N`.")
        new(v, N, similar(v))
    end

    @doc """
        DensityMatrix(v::AbstractVector{<:Real})

    Density matrix representation given Pauli-vector dotted with with vector of coefficients
    """
    function DensityMatrix(v::AbstractVector{<:Real})
        N = intlog2(length(v)) ÷ 2
        length(v) == 4^N || error("Expected length of coefficient vector for density matrix is `4^N` for some integer `N`.")
        new(v, N, similar(v))
    end
end


"""
    matrix(ρ::DensityMatrix)

Matrix representation of density matrix `ρ`.
"""
function matrix(ρ::DensityMatrix)
    # Pauli matrix basis (including identity matrix)
    halfpauli = SparseMatrixCSC{ComplexQ,Int}[
        0.5 * sparse(one(ComplexQ)*I, 2, 2),
        0.5 * sparse_matrix(X),
        0.5 * sparse_matrix(Y),
        0.5 * sparse_matrix(Z),
    ]
    mat = zeros(ComplexQ, 2^ρ.N, 2^ρ.N)
    m = quaternary_digits(ρ.N, 0)
    for i in 1:4^ρ.N
        mat += ρ.v[i] * kron([halfpauli[p + 1] for p in reverse(quaternary_digits!(m, i - 1))]...)
    end
    return mat
end


"""
    density_from_statevector(ψ::Vector)

Construct density matrix |ψ⟩⟨ψ| (Pauli representation) corresponding to quantum state `ψ`.
```jldoctest
julia> ψ = 1/sqrt(2)[1, -1]; # create ket - state.
julia> ρ = density_from_statevector(ψ)
DensityMatrix([0.9999999999999998, -0.9999999999999998, 0.0, 0.0], 1)
```
"""
function density_from_statevector(ψ::Vector)
    N = intlog2(length(ψ))
    2^N == length(ψ) || error("Length of input vector must be a power to 2.")
    pauli = SparseMatrixCSC{ComplexQ,Int}[
        sparse(one(ComplexQ)*I, 2, 2),
        sparse_matrix(X),
        sparse_matrix(Y),
        sparse_matrix(Z),
    ]
    v = zeros(4^N)
    m = quaternary_digits(N, 0)
    for i in 1:4^N
        v[i] = real(dot(ψ, kron([pauli[p + 1] for p in reverse(quaternary_digits!(m, i - 1))]...) * ψ))
    end
    return DensityMatrix(v, N)
end


"""
    density_from_matrix(ρmat::AbstractMatrix)

Construct density matrix in Pauli representation from matrix `ρmat`.
```jldoctest
julia> ρmat = 0.15 .*[ 1 -1; -1  1] + 0.35 .* [1 1; 1 1]; # create mixed state of + and -.
julia> ρ = density_from_matrix(ρmat)
DensityMatrix([1.0, 0.39999999999999997, 0.0, 0.0], 1)
```
"""
function density_from_matrix(ρmat::AbstractMatrix)
    N = intlog2(size(ρmat, 1))
    (2^N == size(ρmat, 1) && 2^N == size(ρmat, 2)) || error("Input must be a square `2^N × 2^N` matrix.")
    pauli = SparseMatrixCSC{ComplexQ,Int}[
        sparse(one(ComplexQ)*I, 2, 2),
        sparse_matrix(X),
        sparse_matrix(Y),
        sparse_matrix(Z),
    ]
    v = zeros(4^N)
    m = quaternary_digits(N, 0)
    for i in 1:4^N
        v[i] = real(tr(kron([pauli[p + 1] for p in reverse(quaternary_digits!(m, i - 1))]...) * ρmat))
    end
    return DensityMatrix(v, N)
end


"""
    pauli_group_matrix(ipauli)

Construct a Pauli group matrix (N-fold tensor products of identity and Pauli matrices)
encoded as `DensityMatrix`. `ipauli[j]` specifies the j-th matrix as integer 0,...,3.
"""
function pauli_group_matrix(ipauli::AbstractVector{<:Integer})
    i4 = Matrix{Int}(I, 4, 4)
    DensityMatrix(kron([2 * i4[i + 1, :] for i in reverse(ipauli)]...), length(ipauli))
end

function pauli_group_matrix(paulistring::String)
    pd = Dict('I' => 0, 'X' => 1, 'Y' => 2, 'Z' => 3)
    pauli_group_matrix([pd[s] for s in paulistring])
end
