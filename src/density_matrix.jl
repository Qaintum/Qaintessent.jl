
"""
    DensityMatrix{N}

Density matrix, represented with respect to identity and Pauli basis by storing
4 real coefficients for each qubit, i.e., a real vector of length `4^N` for `N` qubits.
"""
struct DensityMatrix{N}
    "coefficients with respect to identity and Pauli basis"
    v::AbstractVector{<:Real}

    function DensityMatrix{N}(v::AbstractVector{<:Real}) where {N}
        length(v) == 4^N || error("Expected length of coefficient vector for density matrix is `4^N`.")
        new{N}(v)
    end
end


function matrix(ρ::DensityMatrix{N}) where {N}
    # Pauli matrix basis (including identity matrix)
    pauli = [Matrix{Float64}(I, 2, 2), matrix(X), matrix(Y), matrix(Z)]
    mat = zeros(Complex{eltype(ρ.v)}, 2^N, 2^N)
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        mat += ρ.v[i] * kron([pauli[p+1] for p in pt]...)
    end
    return mat
end


function density_from_statevector(ψ::AbstractVector)
    N = convert(Int, log2(length(ψ)))
    @assert 2^N == length(ψ)
    halfpauli = [sparse(0.5I, 2, 2), sparse(0.5*matrix(X)), sparse(0.5*matrix(Y)), sparse(0.5*matrix(Z))]
    v = zeros(4^N)
    for (i, pt) in enumerate((cartesian_tuples(4, N)))
        v[i] = real(dot(ψ, kron([halfpauli[p+1] for p in pt]...) * ψ))
    end
    return DensityMatrix{N}(v)
end


"""
    pauli_group_matrix(ipauli)

Construct a Pauli group matrix (N-fold tensor products of identity and Pauli matrices)
encoded as `DensityMatrix{N}`. `ipauli[j]` specifies the j-th matrix as integer 0,...,3.
"""
function pauli_group_matrix(ipauli::AbstractVector{<:Integer})
    i4 = Matrix{Int}(I, 4, 4)
    DensityMatrix{length(ipauli)}(kron([i4[i+1, :] for i in reverse(ipauli)]...))
end

function pauli_group_matrix(paulistring::String)
    pd = Dict('I' => 0, 'X' => 1, 'Y' => 2, 'Z' => 3)
    pauli_group_matrix([pd[s] for s in paulistring])
end
