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