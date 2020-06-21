using Test
using TestSetExtensions
using Qaintessent
using LightGraphs

function random_graph(Nn, Ne)
    G = Graph(Nn)
    for j in 1:Ne
        n1 = rand(1:Nn-1)
        n2 = rand(n1:Nn)
        add_edge!(G, n1, n2)
    end
    return G
end

# Interaction graph of circuit with only k-neighbors interactions
function local_circuit_graph(N, k)
    G = Graph(N)
    for i in 1:N-1
        for j in 1:min(k-1, N-i)
            LightGraphs.add_edge!(G, i, i+j)
        end
    end
    return G
end


@testset ExtendedTestSet "tree decomposition" begin

    # test for complete graph
    for n in [10, 25, 50]
        G = complete_graph(n)
        tw, _ = tree_decomposition(G)
        @test tw == n-1
    end

    # test tree decomposition for random graph
    for i in 1:5
        Nn = rand(20:50)
        Ne = rand(2*Nn:10*Nn)
        G = random_graph(Nn, Ne)
        tw, tree, bags = tree_decomposition(G)
        @test is_tree_decomposition(G, tree, bags)
    end

    # test treewidth of local circuit graph (heuristic should generate an optimal tree decomposition)
    for n in 2:5
        IG = local_circuit_graph(10, n)
        tw, _ = tree_decomposition(IG)
        @test tw == n-1
    end
end
