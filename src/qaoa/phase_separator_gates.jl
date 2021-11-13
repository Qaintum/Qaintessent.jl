using Qaintessent
using LinearAlgebra
using SparseArrays: sparse
using Memoize

# Simple struct that represents a graph via its edges
struct Graph
    n::Int # number of vertices (indices 1,...,n)
    edges::Set{Set{Int}} # edges, represented as sets of vertices

    function Graph(n::Integer, edges::Vector{Tuple{T, T}}) where T <: Integer
        n >= 1 || throw(DomainError("n must be a positive integer"))

        # Turn into set of sets
        edge_set = Set(Set.(edges))

        # Verify that all edges are valid
        all(edge -> edge ⊆ 1:n && length(edge) == 2, edge_set) || throw(ArgumentError("Some edges have invalid endpoints"))
        new(n, edge_set)
    end
end

"""
    Phase separation gate for Max-κ-colorable subgraph QAOA mapping.
    Represents the objective function which counts the number of invalid
    edges (i.e. between vertices of the same color).
    Implemented the for one-hot encoding.

``U_{P}(\\gamma) = e^{-i \\gamma H_{P}}``
``H_{P} = \\sum_{\\{u, v\\} \\in E} \\sum_{a=1}^{\\kappa} Z_{u, a} Z_{v, a}``

Reference:\n
    Stuart Hadfield, Zhihui Wang, Bryan O'Gorman, Eleanor G. Rieffel, Davide Venturelli and Rupak Biswas\n
    From the Quantum Approximate Optimization Algorithm to a Quantum Alternating Operator Ansatz\n
    Algorithms 12.2 (2019), p.34
"""
struct MaxKColSubgraphPhaseSeparationGate <: AbstractGate 
    # use a reference type (array with 1 entry) for compatibility with Flux
    γ::Vector{Float64} 
    κ::Int # the number of possible colors
    graph::Graph # the underlying graph which should be colored

    function MaxKColSubgraphPhaseSeparationGate(γ::Float64, κ::Integer, graph::Graph)
        κ > 0 || throw(ArgumentError("Parameter `κ` must be a positive integer."))
        length(graph.edges) > 0 || throw(ArgumentError("Graph `graph` must have at least one edge."))
        new([γ], κ, graph)
    end
end

@memoize function max_k_col_subgraph_phase_separation_hamiltonian(graph::Graph, κ::Int)
    z = matrix(Z)
    
    # Implementation of Eq. (17)
    # one-hot encoding: n * κ vector. Index (a-1)*n + (b-1) corresponds to vertex a, color b.
    H_P_enc = sum(
        kron((color == a && vertex ∈ edge ? z : I(2) for vertex ∈ 1:graph.n for color ∈ 1:κ)...) # Z_{u,a} Z_{v,a}
        for a ∈ 1:κ for edge ∈ graph.edges # Σ_{(u,v) = edge ∈ E} Σ_{a=1..κ}
    )
    return Diagonal(H_P_enc)
end

function Qaintessent.matrix(g::MaxKColSubgraphPhaseSeparationGate)
    # Calculate the hamiltonian (Eq. 17)
    H_P_enc = max_k_col_subgraph_phase_separation_hamiltonian(g.graph, g.κ)

    # Implementation of one-hot phase seperator, A.3.1, p. 34
    U_P = exp(-im * g.γ[] * H_P_enc)

    U_P
end

Qaintessent.adjoint(g::MaxKColSubgraphPhaseSeparationGate) = MaxKColSubgraphPhaseSeparationGate(-g.γ[], g.κ, g.graph)

Qaintessent.sparse_matrix(g::MaxKColSubgraphPhaseSeparationGate) = sparse(matrix(g))

# Number of wires (κ * n in the one-hot encoding)
Qaintessent.num_wires(g::MaxKColSubgraphPhaseSeparationGate)::Int = g.κ * g.graph.n
