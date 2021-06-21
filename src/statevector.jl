"""
    Statevector

statevector representation of quantum state
```
julia> s = Statevector(3)
```
"""
struct Statevector
    state::Vector{ComplexF64}
    N::Int
    function Statevector(N::Int)
        N > 0 || error("Statevector must contain at least 1 qubit")
        state = zeros(2^N)
        state[1] = 1
        new(state, N)
    end

    function Statevector(state::Vector{ComplexF64})
        ispow2(length(state)) || error("Statevector must be a length that is a power of two")
        norm(state) == 1 || error("Statevector must have norm 1")
        N = log(2, length(state))
        new(state, N)
    end
end

Base.size(s::Statevector) = s.N