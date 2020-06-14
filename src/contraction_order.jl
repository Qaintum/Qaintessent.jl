using LightGraphs
using PyCall: pyimport

function ⊂(A,B)
    for a in A
        (a ∉ B) ? (return false) : nothing
    end
    return true
end

"""
    simple_dag(cgc)

Simple DAG of a CircuitGateChain using JuliaGraphs
"""
function simple_dag(cgc::CircuitGateChain{N}) where N
    G = Graph()
    node_info = []
    for i in 1:N
        add_vertex!(G)
        push!(node_info, "q$i")
    end
    qubit_node = [1:N...]
    for moment in cgc
        for gate in moment
            add_vertex!(G)
            push!(node_info, gate)
            for i in gate.iwire
                LightGraphs.add_edge!(G, qubit_node[i], nv(G))
                qubit_node[i] = nv(G)
            end
        end
    end
    return G, node_info
end

"""
    line_graph(G)

Line graph of G.
"""
function line_graph(G::Graph)
    LG = Graph()
    node_info = []

    # add nodes
    for i in 1:nv(G)
        neigh = G.fadjlist[i]
        for j in neigh
            # add the node if not present
            edge = (i < j)  ? (i, j) : (j, i)
            if edge ∉ node_info
                add_vertex!(LG)
                push!(node_info, edge)
            end
        end
    end

    # add edges
    for (i, e1) in enumerate(node_info)
        for (j, e2) in enumerate(node_info[i+1:end])
            common = e1 ∩ e2
            if length(common) > 0
                LightGraphs.add_edge!(LG, i, j + i)
            end
        end
    end
    return LG, node_info
end


"""
    interaction_graph(cgc)

Interaction graph of a CircuitGateChain.
"""
function interaction_graph(cgc::CircuitGateChain{N}) where N
    G = Graph(N)
    for moment in cgc
        for cg in moment
            for (j, i1) in enumerate(cg.iwire)
                for i2 in cg.iwire[1:j-1]
                    LightGraphs.add_edge!(G, i1, i2)
                end
            end
        end
    end
    return G
end

"""
    lacking_for_clique_neigh(G, i)

Finds the edges that are lacking to `G` for the neighborhood of vertex `i`
to be a clique. Returns this number of edges and the edges themselves.
"""
function lacking_for_clique_neigh(G, i)
    neigh = G.fadjlist[i]
    lacking = Tuple{Int,Int}[]
    for (j, i1) in enumerate(neigh)
        for i2 in neigh[1:j-1]
            has_edge(G, i1, i2) ? nothing : push!(lacking, (i1, i2))
        end
    end
    return length(lacking), lacking
end

"""
    rem_vertex_fill!(G, i, ordering, vertex_label)

Removes vertex `i` of `G`, pushes it to `ordering` and updates `vertex_label`
(LightGraphs moves the vertices when one is removed).
"""
function rem_vertex_fill!(G, i, ordering, vertex_label)
    push!(ordering, vertex_label[i])
    rem_vertex!(G, i)
    v = pop!(vertex_label)
    if i ≤ nv(G)
        vertex_label[i] = v
    end
    nothing
end


"""
    min_fill_ordering(G)

Find an ordering of the vertices of `G` using the min-fill heuristic
(cfr. [Bodlaender, _Discovering Treewidth_](http://webdoc.sub.gwdg.de/ebook/serien/ah/UU-CS/2005-018.pdf))
"""
function min_fill_ordering(G::Graph)
    H = copy(G)
    ordering = []
    vertex_label = collect(1:nv(H))
    while nv(H) > 0
        success = false
        for i in nv(H):-1:1
            if length(H.fadjlist[i]) == 0
                rem_vertex_fill!(H, i, ordering, vertex_label)
                success = true
            end
        end
        for i in nv(H):-1:1
            if length(H.fadjlist[i]) == 1
                rem_vertex_fill!(H, i, ordering, vertex_label)
                success = true
            end
        end

        if success == false
            degrees = length.(H.fadjlist)
            J = sortperm(degrees)

            found_clique = false
            v = 0
            best_n_lacking = Inf
            best_lacking = Tuple{Int,Int}[]
            for i in 1:nv(H)
                j = J[i]
                n_lacking, lacking = lacking_for_clique_neigh(H, j)
                if n_lacking == 0
                    rem_vertex_fill!(H, j, ordering, vertex_label)
                    found_clique = true
                    break
                elseif n_lacking < best_n_lacking
                    v = j
                    best_n_lacking = n_lacking
                    best_lacking = lacking
                end
            end
            # if running until here, remove v
            if ! found_clique
                for e in best_lacking
                    LightGraphs.add_edge!(H, e[1], e[2])
                end
                rem_vertex_fill!(H, v, ordering, vertex_label)
            end
        end
    end
    ordering
end

"""
    triangulation(G, ordering)

For each vertex add edges between its higher numbered neighbors
"""
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

"""
    tree_decomposition(G)

Finds a tree decomposition of G using the min-fill heuristic. Returns the
width of the decomposition, the tree and the bags.
"""
function tree_decomposition(G::Graph)
    ordering = min_fill_ordering(G)
    H = triangulation(G, ordering)
    up_neighs = [H.fadjlist[ordering[i]] ∩ ordering[i+1:end] for i in 1:nv(H)]
    # first bag is the first node within a maximum clique following the ordering
    clique_neigh_idx = findfirst(length.(up_neighs) .== nv(H)-1:-1:0)
    first_bag = up_neighs[clique_neigh_idx] ∪ ordering[clique_neigh_idx]
    decomp = Graph(1)
    bags = [first_bag]
    tw = length(first_bag) - 1

    for i in clique_neigh_idx-1:-1:1
        neigh = up_neighs[i]

        # find a bag all neighbors are in
        old_bag_idx = 0
        for (j,bag) in enumerate(bags)
            if neigh ⊂ bag
                old_bag_idx = j
                break
            end
        end
        (old_bag_idx == 0) && (old_bag_idx = 1) # no old_bag was found: just connect to the first_bag

        # create new node
        add_vertex!(decomp)
        new_bag = [neigh; ordering[i]]
        push!(bags, new_bag)

        # update treewidth
        tw = max(tw, length(new_bag)-1)

        # add edge to decomposition
        LightGraphs.add_edge!(decomp, old_bag_idx, nv(decomp))

     end
     return tw, decomp, bags
end

"""
    is_tree_decomposition(G, tree, bags)

Checks if `(tree, bags)` forms a tree decomposition of `G`.
"""
function is_tree_decomposition(G, tree, bags)
    # 1. check union property
    println(sort(unique(vcat(bags...))))
    if (unique(sort(vcat(bags...))) != collect(1:nv(G)) )
        @warn ("Union on bags is not equal to union of vertices")
        return false
    end

    # 2. check edge property
    for edge in edges(G)
        edge_found = false
        e = (edge.src, edge.dst)
        for B in bags
            edge_found = edge_found | (e ⊂ B)
        end
        if ! edge_found
            @warn("Some edge not found in any bag")
            return false
        end
    end

    # TODO: 3. check tree property
    return true
end


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
    G_star = G#nx.line_graph(G)
    tw, tree = graph_algorithms_approx.treewidth_min_fill_in(G_star)
    return tw, tree
end

"""
    py_contraction_order(cgc)

Constructs a near optimal contracting order calling `tree_decomposition`, following
Theorem 4.6 of https://arxiv.org/abs/quant-ph/0511069.
"""
function py_contraction_order(cgc::CircuitGateChain{N}) where N
    H, edges = line_graph(simple_dag(cgc)[1])
    tw, tree = py_tree_decomposition(H)
    ordering = []
    while length(tree.nodes) > 1
        bags = collect(tree.nodes)
        ind_leaf = findfirst(tree.degree.(bags) .== 1)
        leaf = bags[ind_leaf]
        neighs = collect(tree.neighbors(leaf)); neigh = neighs[1];
        set_leaf = Set(leaf)
        set_parent = Set(neigh)
        edge_diff = setdiff(set_leaf, set_parent)
        #println(edge_diff)
        for edge in edge_diff
            push!(ordering,edge)
        end
        tree.remove_node(leaf)
    end
    edges[ordering.+1]
end

function contraction_order(cgc::CircuitGateChain{N}) where N
    H, edges = line_graph(simple_dag(cgc)[1])
    tw, (tree, bags) = tree_decomposition(H)
    contr_order = []
    degrees = length.(tree.fadjlist)
    while maximum(degrees) > 0
        leaf_idx = findfirst(degrees .== 1)
        leaf = bags[leaf_idx]
        neigh_idx = tree.fadjlist[leaf_idx][1] # only neighbor of leaf
        neigh = bags[neigh_idx]
        edge_diff_idx = setdiff(leaf, neigh)
        println(edge_diff_idx)
        for i in edge_diff_idx
            push!(contr_order, edges[i])
        end

        # update tree and bags
        rem_vertex!(tree, leaf_idx)
        if leaf_idx ≤ nv(tree)-1 # removed other than last vertex
            moved_bag = pop!(bags)
            bags[leaf_idx] = moved_bag

        end

        degrees = length.(tree.fadjlist)
    end
    contr_order
end
