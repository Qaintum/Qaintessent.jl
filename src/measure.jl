using StatsBase
using Random
using Memoize

"""
    weights(state::Vector{ComplexF64})

extended StatsBase weights to allow complex values
"""
function StatsBase.weights(state::Vector{ComplexF64})
    weights(norm.(state).^2)
end

"""
    measure(state:Vector{ComplexF64}, n::Int)

sample 'n' possible quantum state after measurement.
"""
function measure(state::Vector{ComplexF64}, shots::Int=1)
    ispow2(length(state)) || error("Statevector is not a viable quantum state (length $(length(state) )")
    norm(state) ≈ 1 || error("Statevector is not a viable statevector (norm is $(norm(state))")
    N = Int(log(2, length(state)))

    if shots == 0 || shots == 1
        s = sample(Random.GLOBAL_RNG, [0:(2^N)-1 ... ], weights(state))
        return digits(Int(s), base=2, pad=N) |> reverse
    end

    s = sample(Random.GLOBAL_RNG, [0:(2^N)-1 ... ], weights(state), shots)
    sum(reverse.(digits.(Int.(s), base=2, pad=N))) ./ shots
end