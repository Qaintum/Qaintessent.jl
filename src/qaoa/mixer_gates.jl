using Qaintessent
using LinearAlgebra: I
using SparseArrays: sparse
using Memoize

# Pauli X and Y matrices
x = matrix(X)
y = matrix(Y)

# Utility function for a XY mixer
# Implements X_a X_{a+1} + Y_a Y_{a+1}, or more generally (⊗_{i ∈ xy_indices} X_i) + (⊗_{i ∈ xy_indices} Y_i)
function XY_sum(xy_indices::Vector{Int64}, d::Int64)::Matrix{ComplexF64}
    # passing generator into kron via varargs syntax
    return kron((i ∈ xy_indices ? x : I(2) for i ∈ 1:d)...
        ) + kron((i ∈ xy_indices ? y : I(2) for i ∈ 1:d)...)
end

"""
    r-nearby values single-qudit mixer gate, which acts on a single qudit (but implemented 
    here not for qudits, but the one-hot encoding)

``U_{r\\text{-NV}}(\\beta) = e^{-i \\beta H_{r\\text{-NV}}}``
``H_{r\\text{-NV}} = \\sum_{i=1}^r \\left(\\breve{X}^i + \\left(\\breve{X}^\\dagger\\right)^i\\right)``

Reference:\n
    Stuart Hadfield, Zhihui Wang, Bryan O'Gorman, Eleanor G. Rieffel,Davide Venturelli and Rupak Biswas\n
    From the Quantum Approximate Optimization Algorithm to a Quantum Alternating Operator Ansatz\n
    Algorithms 12.2 (2019), p.34
"""
struct RNearbyValuesMixerGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    β::Vector{Float64}
    r::Integer # this is the r in r-Nearby-values
    d::Integer # d (= κ) = number of colors

    function RNearbyValuesMixerGate(β::Float64, r::Integer, d::Integer)
        (r > 0 && d > 0) || throw(ArgumentError("Parameters r and d must be positive integers."))
        (r <= d-1) || throw(ArgumentError("Parameter d must be between 1 and d-1."))
        new([β], r, d)
    end
end

# Compute the r-NV mixer Hamiltonian
@memoize function r_nearby_values_hamiltonian_onehot(g::RNearbyValuesMixerGate)
    H_rNV = sum(XY_sum([a, ((a + j - 1) % g.d) + 1], g.d) 
        for a in 1:g.d for j in 1:g.r)

    return H_rNV
end

# Compute the r-NV mixer gate matrix: Implementation of Eq. (6)
function Qaintessent.matrix(g::RNearbyValuesMixerGate)
    H_rNV = r_nearby_values_hamiltonian_onehot(g)
    U_rNV = exp(-im * g.β[] * H_rNV)

    return U_rNV
end

# Adjoint of r-NV = r-NV gate with negated β parameter
Qaintessent.adjoint(g::RNearbyValuesMixerGate) = RNearbyValuesMixerGate(-g.β[], g.r, g.d)

Qaintessent.sparse_matrix(g::RNearbyValuesMixerGate) = sparse(matrix(g))

# Number of wires (= g.d since we use the one-hot encoding canonically)
Qaintessent.num_wires(g::RNearbyValuesMixerGate)::Int = g.d

"""
    Parity single-qudit ring mixer gate, which acts on a single qudit (but implemented 
    here not for qudits, but the one-hot encoding)

``U_{\\text{parity}}(\\beta) = U_{\\text{last}}(\\beta) U_{\\text{even}}(\\beta) U_{\\text{odd}}(\\beta)``
``U_{\\text{odd}}(\\beta) = \\prod_{a~\\text{odd}, a \\neq d} e^{-i \\beta (X_a X_{a+1} + Y_a Y_{a+1})}``
``U_{\\text{even}}(\\beta) = \\prod_{a~\\text{even}} e^{-i \\beta (X_a X_{a+1} + Y_a Y_{a+1})}``
``U_{\\text{last}}(\\beta) = e^{-i \\beta (X_d X_1 + Y_d Y_1)} ~\\text{if}~ d ~\\text{is odd,}~ I ~\\text{otherwise.}``

The formulas require some interpretation. Assumptions for this implementation are:
- the indices a start at 1 (not at zero like elsewhere in the paper)
- X_{a+1} and Y_{a+1} actually means X_1 and Y_1 if a = d

Reference:\n
    Stuart Hadfield, Zhihui Wang, Bryan O'Gorman, Eleanor G. Rieffel, Davide Venturelli and Rupak Biswas\n
    From the Quantum Approximate Optimization Algorithm to a Quantum Alternating Operator Ansatz\n
    Algorithms 12.2 (2019), equations (7) - (10), p. 11
"""
struct ParityRingMixerGate <: AbstractGate
    β::Vector{Float64}
    d::Integer # d (= κ) = number of colors
    is_adjoint::Bool # if yes, the order in the product is reversed

    function ParityRingMixerGate(β::Float64, d::Integer; is_adjoint::Bool = false)
        d > 0 || throw(ArgumentError("Parameter d must be a positive integer."))
        new([β], d, is_adjoint)
    end
end

# Compute the parity ring mixer gate matrix: Implements Eq. (8)
function Qaintessent.matrix(g::ParityRingMixerGate)
    # assumption: by a ≠ n, the paper actually means a ≠ d.
    # assumption: by X_a for a = d+1, the paper means X_1.
    U_odd = prod([exp(-im * g.β[] * XY_sum([a, a+1], g.d)) for a ∈ 1:2:(g.d - 1)], init=I)
    U_even = prod([exp(-im * g.β[] * XY_sum([a, a < g.d ? (a+1) : 1], g.d)) for a ∈ 2:2:g.d], init=I)

    # Implements Eq. (9)
    U_last =  isodd(g.d) ? exp(-im * g.β[] * XY_sum([g.d, 1], g.d)) : I

    # Implements Eq. (7) (U_parity)
    if !g.is_adjoint
        return U_last * U_even * U_odd
    else
        return U_odd * U_even * U_last
    end
end

# Adjoint of parity ring mixer -> negated β parameter, reversed order in the product
Qaintessent.adjoint(g::ParityRingMixerGate) = ParityRingMixerGate(-g.β[], g.d, is_adjoint = !g.is_adjoint)

Qaintessent.sparse_matrix(g::ParityRingMixerGate) = sparse(matrix(g))

# Number of wires (= g.d since we use the one-hot encoding canonically)
Qaintessent.num_wires(g::ParityRingMixerGate)::Int = g.d

"""
    Partition single-qudit mixer gate (implemented not for qudits, but the one-hot encoding)

``U_{\\mathcal{P}-r-\\text{NV}}(\\beta) = U_{P_p-\\text{XY}}(\\beta) \\dots U_{P_1-\\text{XY}}(\\beta)``
``U_{P-\\text{XY}}(\\beta) = \\prod_{\\{a, b\\} \\in P} e^{-i \\beta (\\ket{a}\\bra{b} + \\ket{b}\\bra{a})}``

Note: the time evolution term is effectively implemented with an additional factor of two in the exponent:
``e^{-i \\beta 2 (\\ket{a}\\bra{b} + \\ket{b}\\bra{a})}``
to be consistent with the XY gates used for the other mixers.

Reference:\n
    Stuart Hadfield, Zhihui Wang, Bryan O'Gorman, Eleanor G. Rieffel, Davide Venturelli and Rupak Biswas\n
    From the Quantum Approximate Optimization Algorithm to a Quantum Alternating Operator Ansatz\n
    Algorithms 12.2 (2019), equations (11) - (12), p. 12
"""
struct PartitionMixerGate <: AbstractGate
    β::Vector{Float64}
    d::Int64
    partition::Vector{Vector{Tuple{Int, Int}}} # the Tuple stops Flux.@functor from misinterpreting these as params

    function PartitionMixerGate(β::Float64, d::Int64, partition::Vector{Vector{Tuple{Int, Int}}})
        d > 0 || throw(ArgumentError("Parameter d must be a positive integer."))

        # check that no duplicate indices occur in a part (s.t. the XY mixers within one part commute)
        for partition_part ∈ partition
            part_indices = union(reduce(vcat, collect.(partition_part)))
            part_indices ⊆ 1:d || throw("Indices in partition must be between 1 and d.")
            length(part_indices) == 2 * length(partition_part) ||
                throw("No index must occur more than once within each partition part in `partition`.")
        end

        new([β], d, partition)
    end
end

# Compute the partition mixer Hamiltonians
@memoize function partition_mixer_hamiltonians(g::PartitionMixerGate)::Vector{Matrix{ComplexF64}}
    hamiltonians = Matrix{ComplexF64}[]

    # Implements parts of Eqs. (11), (12)
    # iterate through partition
    for partition_part ∈ g.partition
        # we can represent each part by one Hamiltonian, because the individual XY gates commute
        # and therefore exp(-iβ ∑H_{a,b}) = ∏exp(-iβ H_{a,b})
        H_part = sum(XY_sum([a, b], g.d) for (a, b) ∈ partition_part)
        push!(hamiltonians, H_part)
    end

    return hamiltonians
end

# Compute the partition mixer gate matrix, but allowing to pass a different β
# than is stored in the gate struct (useful for backward pass)
function partition_mixer_gate_matrix(g::PartitionMixerGate, β::Float64)
    hamiltonians = partition_mixer_hamiltonians(g)

    # Implements parts of Eqs. (11), (12)
    # reverse terms in product to have matrix application from right to left
    Us = exp.(-im * β * hamiltonians)
    return prod(reverse(Us))
end

Qaintessent.matrix(g::PartitionMixerGate) = partition_mixer_gate_matrix(g, g.β[])

# Adjoint of partition ring mixer -> negated β parameter, reversed order in the product
Qaintessent.adjoint(g::PartitionMixerGate) = PartitionMixerGate(-g.β[], g.d, reverse(g.partition))

Qaintessent.sparse_matrix(g::PartitionMixerGate) = sparse(matrix(g))

# Number of wires (= g.d since we use the one-hot encoding canonically)
Qaintessent.num_wires(g::PartitionMixerGate)::Int = g.d

"""
    Mixer hamiltonian gate for the WSQAOA

``U_{WS}(\\beta) = e^{-i \\beta H_{WS}}``
``H_{WS} =  \\sum_{i=1}^r \\breve{R_{Y}}(\\theta_i)\\breve{R_{Z}}(-2\\beta)\\breve{R_{Y}}(-\\theta_i) where r is the number of qubits``

Reference:\n
    Egger, Daniel J., Jakub Mareček, en Stefan Woerner\n
    “Warm-starting quantum optimization”\n
    Quantum 5 (Junie 2021): 479. https://doi.org/10.22331/q-2021-06-17-479
"""

struct WSQAOAMixerGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    β::Vector{Float64}
    c_opt::Vector{Float64} # solution of the continous relaxation of the initial problem (e.g: MaxCut)
    ε::Float64 # regularization parameter ∈ [0, 0.5], prevents reachability issues
    init_state_randomized::Bool # Is the initial state going to be randomized before calculation?
                   # if yes, then off-diagonal elements of the H_{i} "would be multiplied by -1
                   # to be able to both represent solutions of the random-hyperplane 
                   # rounding as well as deviate from them" (Sec. 2.3, p.4)

    function WSQAOAMixerGate(β::Float64, c_opt::Vector{Float64}, ε::Float64, init_state_randomized::Bool)
        (ε >= 0  && ε <= 0.5) || throw(ArgumentError("Parameters ε be between 0 and 0.5"))
        new([β], c_opt, ε, init_state_randomized)
    end
end

# Compute the WSQAOA mixer Hamiltonian
@memoize function wsqaoa_mixer_hamiltonian(g::WSQAOAMixerGate)::Matrix{ComplexF64}
    H_dim = 2 ^ length(g.c_opt)
    H_wsqaoa = zeros(ComplexF64, H_dim, H_dim)

    # Implements Hamiltonian Part of Eq. (2)
    # iterate through the c_opt (number of qubits)
    for c_i_index in eachindex(g.c_opt)
        # perform regularization
        c_i = g.c_opt[c_i_index]
        c_i = c_i <= g.ε ? g.ε : c_i
        c_i = c_i >= (1 - g.ε) ? (1 - g.ε) : c_i
        H_c_part_off_diagonal_val = -2 * sqrt(c_i * (1 - c_i))
        if g.init_state_randomized
            H_c_part_off_diagonal_val = -H_c_part_off_diagonal_val
        end       
        H_c_part = [(2 * c_i - 1) H_c_part_off_diagonal_val; H_c_part_off_diagonal_val (1 - 2 * c_i)]
        H_i = kron((i == c_i_index ? H_c_part : I(2) for i ∈ 1:length(g.c_opt))...)
        H_wsqaoa += H_i
    end

    return H_wsqaoa
end

# Compute the WSQAOA mixer gate matrix: Implementation of complete Eq. (2)
function Qaintessent.matrix(g::WSQAOAMixerGate)
    U_wsqaoa = [1]

    # Implements mixer gate according to Fig. (2)
    for c_i_index in eachindex(g.c_opt)
        # perform regularization
        c_i = g.c_opt[c_i_index]
        c_i = c_i <= g.ε ? g.ε : c_i
        c_i = c_i >= (1 - g.ε) ? (1 - g.ε) : c_i
        theta_i = 2 * asin(sqrt(c_i))
        β = g.β[]
        # reverse order if randomized
        if g.init_state_randomized
            theta_i = -theta_i
            β = -β
        end
        U_i = matrix(RyGate(theta_i)) * matrix(RzGate(-2 * g.β[])) * matrix(RyGate(-theta_i))
        U_wsqaoa = kron(U_wsqaoa, U_i)
    end

    return U_wsqaoa
end

# Adjoint of WSQAOA mixer -> negated β parameter, reversed order in the product
# here we take advantage of the fact that randomization
# reverses the order of rotations (Sec. 2.3, p.4)
Qaintessent.adjoint(g::WSQAOAMixerGate) = WSQAOAMixerGate(-g.β[], g.c_opt, g.ε, g.init_state_randomized)

Qaintessent.sparse_matrix(g::WSQAOAMixerGate) = sparse(matrix(g))

# Number of wires (= g.d since we use the one-hot encoding canonically)
Qaintessent.num_wires(g::WSQAOAMixerGate)::Int = length(g.c_opt)

# Classical Rx Mixer Gate typically used for MaxCut, implemented for comparison
struct RxMixerGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    β::Vector{Float64}
    n::Int

    function RxMixerGate(β::Float64, n::Int)
        new([β], n)
    end
end

@memoize function rx_mixer_hamiltonian(g::RxMixerGate)::Matrix{ComplexF64}
    H_dim = 2 ^ g.n
    H_rx = zeros(ComplexF64, H_dim, H_dim)
    for i in 1:g.n
        H_rx_i = kron((i == j ? matrix(X) : I(2) for j ∈ 1:g.n)...)
        H_rx += -H_rx_i
    end

    return H_rx
end

function Qaintessent.matrix(g::RxMixerGate)
    U_rx = [1.0]

    for i in 1:g.n
        U_i = matrix(RxGate(2.0 * g.β[]))
        U_rx = kron(U_rx, U_i)
    end
    
    return U_rx
end

Qaintessent.adjoint(g::RxMixerGate) = RxMixerGate(-g.β[], g.n)

Qaintessent.sparse_matrix(g::RxMixerGate) = sparse(matrix(g))

Qaintessent.num_wires(g::RxMixerGate)::Int = g.n

