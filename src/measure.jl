using StatsBase
using Random
using Memoize

"""
    weights(state::Vector{ComplexQ})

extended StatsBase weights to allow complex values
"""
function StatsBase.weights(state::Vector{ComplexQ})
    weights(norm.(state).^2)
end

"""
    measure(state:Vector{ComplexQ}, n::Int)

sample 'n' possible quantum state after measurement.
"""
function measure(state::Vector{ComplexQ}, shots::Int=1)

    (shots > 0) || error("Shots has to be a natural integer above 0.")	
    ispow2(length(state)) || error("Statevector is not a viable quantum state (length $(length(state) )")
    norm(state) ≈ 1 || error("Statevector is not a viable statevector (norm is $(norm(state))")
    N = Int(log(2, length(state)))

    # generate samples and count number of measurements of each basis state 
    s = sample(Random.GLOBAL_RNG, [0:(2^N)-1 ... ], weights(state), shots)     
    u=unique(s)
    d=Dict([(i,count(x->x==i,s)) for i in u]) 
    
    return d
end
