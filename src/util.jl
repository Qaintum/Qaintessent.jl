
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
pauli_vector(x, y, z) = [z x - im * y; x + im * y -z]
