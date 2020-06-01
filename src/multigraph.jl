mutable struct MultiGraph
    nodes::AbstractVector{Any}                  # labels of the nodes
    edges::AbstractVector{Tuple{Int, Int, Any}} # edges and their labels
    adj_node_node::AbstractVector{Set{Int}}     # adjacency between nodes
    adj_node_edge::AbstractVector{Set{Tuple{Int, Int, Any}}} # adj. between nodes and edges

    function MultiGraph()
        return new(Any[], Tuple{Int, Int, Any}[], Set{Int}[], Set{Tuple{Int, Int, Any}}[])
    end

    function MultiGraph(nodes, edges)
        G = MultiGraph()
        for n in nodes
            add_node!(G, n)
        end
        for e in edges
            add_edge!(G, e)
        end
        return G
    end
end

Base.length(G::MultiGraph) = length(G.nodes)

function add_node!(G, label)
    push!(G.nodes, label)
    push!(G.adj_node_node,Set{Int}())
    push!(G.adj_node_edge,Set{Tuple{Int, Int, Any}}())
    nothing
end

add_node!(G) = add_node!(G, length(G)+1)

get_node_idx(G, key) = findfirst(broadcast(x -> (x == key), G.nodes))

function add_edge!(G, i, j, key; warn_repeating = true)
    i < j ? nothing : (i,j) = (j,i)
    e = (i, j, key)
    if (e ∈ G.edges)
        warn_repeating && @warn "Edge ($i, $j) with key $(key) already exists."
    else
        push!(G.edges, e)
        push!(G.adj_node_node[i], j)
        push!(G.adj_node_node[j], i)
        push!(G.adj_node_edge[i], e)
        push!(G.adj_node_edge[j], e)
    end
    nothing
end

add_edge!(G, i, j) = add_edge!(G, i, j, 0)

add_edge!(G, tup::Tuple) = add_edge!(G, tup...)

function simple_dag(cgc::CircuitGateChain{N}) where N
    G = MultiGraph()
    for i in 1:N
        add_node!(G, "q$i")
    end
    qubit_node = [1:N...]
    for moment in cgc
        for gate in moment
            add_node!(G, gate)
            for i in gate.iwire
                add_edge!(G, qubit_node[i], length(G), i)
                qubit_node[i] = length(G)
            end
        end
    end
    return G
end

"""
    line_graph(G::MultiGraph)

Creates the line graph of G.
"""
function line_graph(G::MultiGraph)
    LG = MultiGraph()
    for i in 1:length(G)
        node_edges = collect(G.adj_node_edge[i]) # transform to array so that enumerate works always the same
        for (j, e1) in enumerate(node_edges)
            # add the node if not present
            if e1 ∉ LG.nodes
                add_node!(LG, e1)
            end

            # create edges between all co-incident nodes
            for (k, e2) in enumerate(node_edges[1:j-1])
                add_edge!(LG, get_node_idx(LG, e2), get_node_idx(LG, e1), (e2, e1, i), warn_repeating=false)
            end
        end
    end
    return LG
end
