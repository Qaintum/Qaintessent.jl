using PyCall: pyimport

"""
    py_tree_decomposition(adj)

Constructs a near-optimal tree decomposition of the graph whose adjacency matrix is adj.
Calls library `Networkx`.
"""
function py_tree_decomposition(G0::Graph)
    nx = pyimport("networkx")
    graph_algorithms_approx = pyimport("networkx.algorithms.approximation")
    adj = Array(adjacency_matrix(G0))
    G = nx.from_numpy_matrix(adj)
    tw, tree = graph_algorithms_approx.treewidth_min_fill_in(G)
    return tw, tree
end

function py_gplot(PyGraph)
    nx = pyimport("networkx")
    plt = pyimport("matplotlib.pyplot")
    plt.figure(figsize = (6,9))
    nx.draw_networkx(PyGraph)
    xlim = plt.gca().get_xlim()
    plt.gca().set_xlim(xlim[1] - 0.3, xlim[2]+0.3)
    plt.savefig("graphplot.png",dpi = 300)
end


"""
    py_contraction_order(cgc)

Constructs a near optimal contracting order calling `tree_decomposition`, following
Theorem 4.6 of https://arxiv.org/abs/quant-ph/0511069.
"""
function py_contraction_order(cgc::CircuitGateChain{N}) where N
    H, edges = line_graph(simple_dag(cgc)[1])
    tw, tree = py_tree_decomposition(H)
    contr_order = []
    bags = collect(tree.nodes)
    while length(tree.nodes) > 1
        #println(tree.degree.(bags))
        leaves_idx = findall(tree.degree.(bags) .<= 1)
        leaf_idx = leaves_idx[argmin(length.(bags[leaves_idx]))]
        #ind_leaf = findfirst(tree.degree.(bags) .== 1)
        leaf = bags[leaf_idx]
        neighs = collect(tree.neighbors(leaf)); neigh = neighs[1];
        set_leaf = Set(leaf)
        set_parent = Set(neigh)
        edge_diff = setdiff(set_leaf, set_parent)
        #println(edge_diff .+ 0)
        for edge in edge_diff
            push!(contr_order,edge)
        end
        tree.remove_node(leaf)
        bags = collect(tree.nodes)
    end
    # check that only one bag left
    @assert length(bags) == 1
    for b in bags
        for edge in b
            push!(contr_order,edge)
        end
    end
    edges[contr_order.+1]
end
