const Edge = Tuple{Int64, Int64, Any}

mutable struct MultiGraph
    nodes::AbstractVector{Any}                  # labels of the nodes
    edges::AbstractVector{Edge} # edges and their labels
    adj_node_node::AbstractVector{Set{Int}}     # adjacency between nodes
    adj_node_edge::AbstractVector{Set{Edge}} # adj. between nodes and edges

    function MultiGraph()
        return new(Any[], Edge[], Set{Int}[], Set{Edge}[])
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
    push!(G.adj_node_edge,Set{Edge}())
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

function delete_repeated_edges(G::MultiGraph)
    G = MultiGraph(G.nodes, G.edges)
    total_repeated = Set{Edge}()
    for n in 1:length(G)
        c = 1
        node_edges = collect(G.adj_node_edge[n])
        while c < length(node_edges)
            e1 = node_edges[c]
            repeated = []
            for e2 in node_edges[c+1:end]
                if (e1[1] == e2[1]) & (e1[2] == e2[2])
                    push!(repeated, e2)
                end
            end
            mask = [e ∉ repeated for e in node_edges]
            node_edges = node_edges[mask]
            if length(repeated) > 0
                push!(total_repeated, repeated...)
            end
            c += 1
        end
        G.adj_node_edge[n] = Set{Edge}(node_edges)
    end
    mask = [e ∉ total_repeated for e in collect(G.edges)]
    G.edges = G.edges[mask]
    return G
end



"""
    line_graph(G; ignore_repeated = true)

Creates the line graph of G. If `ignore_repeated` first deletes all
repeated edges of G.
"""
function line_graph(G::MultiGraph, ignore_repeated = true)
    LG = MultiGraph()

    if ignore_repeated
        G = delete_repeated_edges(G)
    end

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
