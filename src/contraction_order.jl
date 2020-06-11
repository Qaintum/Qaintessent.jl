using LightGraphs
using PyCall

function ⊂(A,B)
    for a in A
        (a ∉ B) ? (return false) : nothing
    end
    return true
end

# compute interaction graph
function interaction_graph(cgc::CircuitGateChain{N}) where N
    G = Graph(N)
    for moment in cgc
        for cg in moment
            for (j, i1) in enumerate(cg.iwire)
                for i2 in cg.iwire[1:j-1]
                    LightGraphs.add_edge!(G,i1,i2)
                end
            end
        end
    end
    return G
end

function min_fill_ordering(G::Graph)
    H = copy(G)
    ordering = []
    vertex_label = collect(1:nv(H))
    """graph_plots = []
    locs_x_fixed = rand(N)
    locs_y_fixed = rand(N)"""
    while length(vertices(H)) > 0
        """push!(vertex_label_cum, vertex_label)
        push!(graph_plots, gplot(H,nodelabel = vertex_label,layout = my_layout))"""
        # TODO: check for isolated vertices and leaves before performing the exploration
        v = 0
        best_n_lacking = Inf
        best_lacking = Tuple{Int,Int}[]
        found_leave = false

        for i in vertices(H)
            n_lacking, lacking = lacking_for_clique_neigh(H, i)
            if n_lacking == 0
                push!(ordering, vertex_label[i])
                rem_vertex!(H, i)
                vertex_label = vertex_label[1:length(vertex_label) .!= i]
                found_leave = true
                break
            elseif n_lacking < best_n_lacking
                v = i
                best_n_lacking = n_lacking
                best_lacking = lacking
            end
        end

        # if running until here, remove v
        if ! found_leave
            for e in best_lacking
                LightGraphs.add_edge!(H, e[1], e[2])
            end
            push!(ordering, vertex_label[v])
            rem_vertex!(H, v)
            vertex_label = vertex_label[1:length(vertex_label) .!= v]

        end
        #println(v)
        #println(best_lacking)

    end
    ordering
end

function triangulation(G, ordering)
    H = copy(G)
    for (i, v) in enumerate(ordering)
        neigh = H.fadjlist[v]
        high_neigh = neigh ∩ ordering[i+1:end] # get higher numbered neighbors
        for (j, i1) in enumerate(high_neigh)
            for i2 in high_neigh[1:j-1]
                LightGraphs.add_edge!(H, i1, i2)
            end
        end
    end
    return H
end


function tree_decomposition(G::Graph)
    ordering = min_fill_ordering(G)
    H = triangulation(G, ordering)

    bags = []
    for (i,v) in enumerate(ordering)
        bag = (H.fadjlist[v] ∩ ordering[i+1:end]) ∪ v
        if (length(bag) > 0) & prod([!(bag ⊂ b) for b in bags]) > 0
            push!(bags,bag)
        end
    end
    maximum(length.(bags))-1, bags
end

tree_decomposition(cgc::CircuitGateChain) = tree_decomposition(interaction_graph(cgc)...)


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
    G_star = nx.line_graph(G)
    tw, tree = graph_algorithms_approx.treewidth_min_fill_in(G_star)
    return tw, tree
end


"""
    contraction_order(cgc)

Constructs a near optimal contracting order calling `tree_decomposition`, following
Theorem 4.6 of https://arxiv.org/abs/quant-ph/0511069.
"""
function py_contraction_order(cgc::CircuitGateChain{N}) where N
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
