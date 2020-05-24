using PyCall

"""
    simple_dag(cgc)

Returns the lists nodes and edges of a CircuitGateChain Dag
"""
function simple_dag(cgc::CircuitGateChain{N}) where N
    nodes = "q" .* string.(1:N)
    qubit_node = [1:N...]
    edges = []
    for moment in cgc
        for gate in moment
            push!(nodes,string(gate))
            for i in gate.iwire
                push!(edges,(qubit_node[i], length(nodes)))
                qubit_node[i] = length(nodes)
            end
        end
    end
    return nodes, edges
end

""" Adjacency matrix of a graph"""
function adj_matrix(nodes, edges)
    N_nodes = length(nodes)
    adj = zeros(Int, N_nodes, N_nodes)
    for (i,j) in edges
        adj[i,j] = 1
    end
    return adj + adj'
end

"""
    tree_decomposition(adj)

Constructs a near-optimal tree decomposition of the graph whose adjacency matrix is adj.
Calls library `Networkx`.
"""
function tree_decomposition(adj::AbstractMatrix)
    nx = pyimport("networkx")
    graph_algorithms_approx = pyimport("networkx.algorithms.approximation")
    G = nx.from_numpy_matrix(adj)
    G_star = nx.line_graph(G)
    tw, tree = graph_algorithms_approx.treewidth_min_fill_in(G_star)
    return tw, tree
end

tree_decomposition(cgc::CircuitGateChain) = tree_decomposition(adj_matrix(simple_dag(cgc)...))

"""
    contraction_order(cgc)

Constructs a near optimal contracting order calling `tree_decomposition`, following
Theorem 4.6 of https://arxiv.org/abs/quant-ph/0511069.
"""
function contraction_order(cgc::CircuitGateChain{N}) where N
    tw, tree = tree_decomposition(cgc)
    ordering = []
    while length(tree.nodes) > 1
        bags = collect(tree.nodes)
        ind_leaf = findfirst(tree.degree.(bags) .== 1)
        leaf = bags[ind_leaf]
        neighs = collect(tree.neighbors(leaf)); neigh = neighs[1];
        set_leaf = Set(leaf)
        set_parent = Set(neigh)
        edge_diff = setdiff(set_leaf, set_parent)
        for edge in edge_diff
            push!(ordering,edge)
        end
        tree.remove_node(leaf)
    end
    ordering
end
