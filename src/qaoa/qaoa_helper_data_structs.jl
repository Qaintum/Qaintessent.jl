# Struct that represents a graph via its weighted edges
struct EdgeWeightedGraph
    n::Int # number of vertices (indices 1,...,n)
    edges::Set{Tuple{Set{Int}, Float64}} # edges, represented as tuples of vertices and weights

    function EdgeWeightedGraph(n::Integer, edges::Vector{Tuple{Int, Int, Float64}})
        n >= 1 || throw(DomainError("n must be a positive integer"))
        edge_set = Set([((Set((v1, v2)), w) for (v1, v2, w) in edges)...])
        all((e, w)::Tuple -> e ⊆ 1:n, edge_set) || throw(ArgumentError("Some edges have invalid endpoints"))
        new(n, edge_set)
    end
end

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

# Trivially convert a unweighted graph to a weighted one with edge weight 1
function to_edge_weighted_graph(g::Graph)::EdgeWeightedGraph
    weighted_edge_set = [((v1, v2, 1.0) for (v1, v2) in g.edges)...]
    return EdgeWeightedGraph(g.n, weighted_edge_set)
end

# Calculates the adjaceny
function adjacency_matrix(g::EdgeWeightedGraph)::Matrix{Int}
    adj_mat = zeros(g.n, g.n)
    for ((v1, v2), w) in g.edges
        adj_mat[v1, v2] = w
        adj_mat[v2, v1] = w
    end
    return adj_mat
end
