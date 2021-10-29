using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using Random
using Qaintessent.QAOAHelperDataStructs
using Qaintessent.MaxKColSubgraphQAOA

graphs = [
    Graph(3, [(1,2), (2,3)]),
    Graph(5, [(1,2), (2,3), (3,4), (4,5), (5,1), (1,3)]),
    Graph(4, [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
]

# Utility function to count the properly colored edges in a graph coloring (i.e. endpoints have different colors)
function properly_colored_edges(graph::Graph, coloring::Vector{Int})
    length(coloring) == graph.n || throw(ArgumentError("Length of coloring must equal the number of vertices."))
    return count(coloring[a] != coloring[b] for (a, b) ∈ graph.edges)
end

# Create a state ψ that corresponds to a given coloring
function ψ_from_coloring(n::Int, κ::Int, colors::Vector{Int})::Vector{ComplexF64}
    (n > 0 && κ > 0) || throw(DomainError("Parameters n and κ must be positive integers."))
    colors ⊆ 1:κ || throw(ArgumentError("Parameter `colors` may only contain colors in the range 1:$(κ)."))

    # Create ψ with |1> entries in the indices corresponding to the given colors, |0> elsewhere
    ψ = kron((color == colors[vertex] ? [0.0im, 1] : [1, 0.0im] for vertex ∈ 1:n for color ∈ 1:κ)...)
    ψ
end

# Create a state ψ that corresponds to a given partitioning
function ψ_from_partitioning(n::Int, partitioning::Vector{Int})::Vector{ComplexF64}
    (n > 0) || throw(DomainError("Parameter n must be positive integers."))
    partitioning ⊆ 0:1 || throw(ArgumentError("Parameter `partitioning` may only contain values in the range 0:1."))

    # Create ψ with |1> entries corresponding to the given partitioning
    ψ = kron((1 == partitioning[vertex] ? [0.0im, 1] : [1, 0.0im] for vertex ∈ 1:n)...)
    ψ
end

# Utility function to compute the probabilities of the outcomes represented by a wavefunction ψ, sorted by descending probability
function wavefunction_distribution(ψ::Vector{ComplexF64}; as_bitstrings::Bool = true,
        include_zero = false)::Union{Vector{Tuple{Int, Float64}}, Vector{Tuple{Vector{Int}, Float64}}}
    distribution = [(i-1, abs(amplitude)^2) for (i, amplitude) ∈ enumerate(ψ) if abs(amplitude) > 0 || include_zero]
    if as_bitstrings
        N = Int(log2(length(ψ)))
        distribution = [(digits(i, base=2, pad=N) |> reverse, p) for (i, p) ∈ distribution]
    end

    return sort(distribution, by=(t -> -t[2]))
end

@testset ExtendedTestSet "max-k-col. subgraph QAOA Hamiltonians and gate matrices" begin
    @testset ExtendedTestSet "max-k-col. subgraph QAOA Hamiltonians" begin
        # Test that the phase separation Hamiltonian used in MaxKColSubgraphPhaseSeparationGate
        # implements the objective function correctly.
        @testset "Hamiltonian MaxKColSubgraphPhaseSeparationGate" begin
            κs = [4, 2, 3]

            hamiltonians = max_k_col_subgraph_phase_separation_hamiltonian.(graphs, κs)
            colorings = [rand(1:κ, graph.n) for (κ, graph) ∈ zip(κs, graphs)]

            for (H, graph, κ, coloring) ∈ zip(hamiltonians, graphs, κs, colorings)
                f_val = properly_colored_edges(graph, coloring)
                ψ_coloring = ψ_from_coloring(graph.n, κ, coloring)
                scaling = κ * length(graph.edges) - 4 * f_val # ψ_coloring should be scaled by this factor
                @test H * ψ_coloring ≈ ψ_coloring * scaling
            end
        end

        # Utility functions to perform circular shift on Vectors (that represent bitstrings)
        left_circ_shift(l::Vector) = [l[2:length(l)]; l[1:1]]
        right_circ_shift(l::Vector) = [l[length(l):length(l)]; l[1:length(l)-1]]

        # Shorthand for repeated function application
        # Taken from https://stackoverflow.com/questions/39895672/apply-function-repeatedly-a-specific-number-of-times
        (repeatApply)(f::Function, i::Int) = i==1 ? f : x->(repeatApply(f, i-1))(f(x))

        # Test that the Hamiltonian which the RNearbyValuesMixerGate is based on
        # maps coloring states to a superposition of their neighboring states
        @testset "Hamiltonian RNearbyValuesMixerGate" begin
            ds = [2, 5, 9, 10] # different numbers of colors
            rs = [1, 2, 5]
            nums_samples = [2, 5, 5, 5]

            for (d, num_samples) ∈ zip(ds, nums_samples)
                for r ∈ rs[rs .< d] # only test valid rs
                    # compute the Hamiltonian
                    gate = RNearbyValuesMixerGate(0., r, d)
                    hamiltonian = r_nearby_values_hamiltonian_onehot(gate)

                    # for each color (i.e. |1000..00>, |0100..00>, ..., |0000..01>):
                    for color ∈ randperm(d)[1:num_samples]
                        # compute superposition state bitstrings after applying Hamiltonian
                        ψ = ψ_from_coloring(1, d, [color])
                        ψ_bits = wavefunction_distribution(ψ)[1][1]
                        out_vec = hamiltonian * ψ # not a wavefunction (we apply a Hamiltonian, not a gate)
                        superposition_states = (wavefunction_distribution(out_vec)
                            |> (distr -> filter(d -> !(d[2] ≈ 0), distr))
                            .|> (d -> d[1]))

                        # compute bitstrings of nearby colors by shifting
                        near_values = [
                            [repeatApply(left_circ_shift, k)(ψ_bits) for k ∈ 1:r];
                            [repeatApply(right_circ_shift, k)(ψ_bits) for k ∈ 1:r]
                        ] # repeatApply is repeated function application defined above

                        # r-NV Mixer should map state to superposition of its neighbor colors
                        @test Set(superposition_states) == Set(near_values)
                    end
                end
            end
        end
    end

    # Utility function that returns the Hamming weights of all basis states of ψ
    hamming_weights(ψ) = (wavefunction_distribution(ψ)
    |> (distr -> filter(d -> !(d[2] ≈ 0), distr))
    .|> (d -> sum(d[1])))


    @testset ExtendedTestSet "max-k-col. subgraph QAOA mixer Hamming weights" begin
        # Test that the RNearbyValuesMixerGate preserves the Hamming weight of coloring
        # states (i.e. maps states of Hamming weight one to a superposition of states
        # with Hamming weight one)
        @testset "Hamming weight RNearbyValuesMixerGate" begin
            ds = [2, 5, 9, 10] # different numbers of colors
            rs = [1, 2, 5]
            nums_samples = [2, 5, 5, 5]

            for (d, num_samples) ∈ zip(ds, nums_samples)
                for r ∈ rs[rs .< d] # only test valid rs
                    gate = RNearbyValuesMixerGate(rand() * 10 + 0.1, r, d)
                    U = Qaintessent.matrix(gate)

                    for color ∈ randperm(d)[1:num_samples]
                        ψ = ψ_from_coloring(1, d, [color])
                        ψ_out = U * ψ

                        # All basis states should have Hamming weight 1
                        out_hamming_weights = hamming_weights(ψ_out)
                        @test hamming_weights(ψ_out) == ones(Int64, length(out_hamming_weights))
                    end
                end
            end
        end

        # Test that the ParityRingMixerGate preserves the Hamming weight of coloring
        # states (i.e. maps states of Hamming weight one to a superposition of states
        # with Hamming weight one), if d == 2
        @testset "Hamming weight ParityRingMixerGate (d = 2)" begin
            ds = [2] # different numbers of colors

            for d ∈ ds
                gate = ParityRingMixerGate(rand() * 10 + 0.1, d)
                U = Qaintessent.matrix(gate)

                for color ∈ 1:d
                    ψ = ψ_from_coloring(1, d, [color])
                    ψ_out = U * ψ

                    # All basis states should have Hamming weight 1
                    out_hamming_weights = hamming_weights(ψ_out)
                    @test hamming_weights(ψ_out) == ones(Int64, length(out_hamming_weights))
                end
            end
        end

        # Test that the ParityRingMixerGate preserves the Hamming weight of coloring
        # states (i.e. maps states of Hamming weight one to a superposition of states
        # with Hamming weight one), if d ≥ 3
        @testset "Hamming weight ParityRingMixerGate (d ≥ 3)" begin
            ds = [5, 8, 9] # different numbers of colors
            num_samples = 5

            for d ∈ ds
                gate = ParityRingMixerGate(rand() * 10 + 0.1, d)
                U = Qaintessent.matrix(gate)

                for color ∈ randperm(d)[1:num_samples]
                    ψ = ψ_from_coloring(1, d, [color])
                    ψ_out = U * ψ

                    # All basis states should have Hamming weight 1
                    out_hamming_weights = hamming_weights(ψ_out)
                    @test hamming_weights(ψ_out) == ones(Int64, length(out_hamming_weights))
                end
            end
        end

        # Test that the PartitionMixerGate preserves the Hamming weight of coloring
        # states (i.e. maps states of Hamming weight one to a superposition of states
        # with Hamming weight one)
        @testset "Hamming weight PartitionMixerGate" begin
            ds = [5, 8, 9] # different numbers of colors
            num_samples = 5
            nums_partition_layers = [3, 2, 7]

            for (d, num_partition_layers) ∈ zip(ds, nums_partition_layers)
                # create num_partition_layers partition parts, each of
                # random length and with random entries
                partition = [
                    # create h x 2 matrix with random height 2 from a random permutation of the indices...
                    (reshape(randperm(d)[1:(2*rand(1:div(d, 2)))], (:, 2))
                    # ... then map the rows to tuples
                        |> eachrow .|> Tuple)
                    for _ ∈ 1:num_partition_layers
                ]

                gate = PartitionMixerGate(rand() * 10 + 0.1, d, partition)
                U = Qaintessent.matrix(gate)

                for color ∈ randperm(d)[1:num_samples]
                    ψ = ψ_from_coloring(1, d, [color])
                    ψ_out = U * ψ

                    # All basis states should have Hamming weight 1
                    out_hamming_weights = hamming_weights(ψ_out)
                    @test hamming_weights(ψ_out) == ones(Int64, length(out_hamming_weights))
                end
            end
        end
    end

    # Test that the parity ring mixer is a special case of the partition mixer
    @testset "PartitionMixerGate generalizes ParityRingMixerGate" begin
        ds = [2, 5, 9, 10]

        for d ∈ ds
            β = rand() * 10 + 0.1

            # construct the parity mixer
            parity_mixer = ParityRingMixerGate(β, d)

            # construct the equivalent partition mixer
            odd_part = [(a, a+1) for a ∈ 1:2:(d-1)]
            even_part = [(a, (a % d) + 1) for a ∈ 2:2:d]
            last_part = isodd(d) ? [(d, 1)] : nothing
            partition = isodd(d) ? [odd_part, even_part, last_part] : [odd_part, even_part]
            partition_mixer = PartitionMixerGate(β, d, partition)

            @test Qaintessent.matrix(parity_mixer) ≈ Qaintessent.matrix(partition_mixer)
        end
    end

    @testset ExtendedTestSet "QAOA gates adjoints" begin
        # Test that `MaxKColSubgraphPhaseSeparationGate` has the correct adjoint.
        @testset "adjoint MaxKColSubgraphPhaseSeparationGate" begin
            γs = rand(length(graphs)) * 2π
            κs = [3, 2, 2]

            gates = MaxKColSubgraphPhaseSeparationGate.(γs, κs, graphs)
            for i ∈ 1:length(gates)
                @test Qaintessent.matrix(gates[i]) * Qaintessent.matrix(Base.adjoint(gates[i])) ≈ I
            end
        end

        # Test that `ParityRingMixerGate` has the correct adjoint.
        @testset "adjoint ParityRingMixerGate" begin
            N = 3 # number of test cases
            βs = rand(N) * 2π
            ds = rand(2:9, N)

            gates = ParityRingMixerGate.(βs, ds)
            for i ∈ 1:length(gates)
                @test Qaintessent.matrix(gates[i]) * Qaintessent.matrix(Base.adjoint(gates[i])) ≈ I
            end
        end

        # Test that `RNearbyValuesMixerGate` has the correct adjoint.
        @testset "adjoint RNearbyValuesMixerGate" begin
            βs = rand(4) * 2π
            ds = [2, 3, 6, 9]
            rs = [1, 2, 1, 5]

            gates = RNearbyValuesMixerGate.(βs, rs, ds)
            for i ∈ 1:length(gates)
                @test Qaintessent.matrix(gates[i]) * Qaintessent.matrix(Base.adjoint(gates[i])) ≈ I
            end
        end

        # Test that `PartitionMixerGate` has the correct adjoint.
        @testset "adjoint PartitionMixerGate" begin
            ds = [2, 3, 6, 9]
            nums_partition_layers = [2, 3, 4, 5]

            for (d, num_partition_layers) ∈ zip(ds, nums_partition_layers)
                # create num_partition_layers partition parts, each of
                # random length and with random entries
                partition = [
                    # create h x 2 matrix with random height 2 from a random permutation of the indices...
                    (reshape(randperm(d)[1:(2*rand(1:div(d, 2)))], (:, 2))
                    # ... then map the rows to tuples
                        |> eachrow .|> Tuple)
                    for _ ∈ 1:num_partition_layers
                ]

                β = rand() * 10 + 0.1
                gate = PartitionMixerGate(β, d, partition)
                @test Qaintessent.matrix(gate) * Qaintessent.matrix(Base.adjoint(gate)) ≈ I
            end
        end
    end

    @testset ExtendedTestSet "max-k-col. subgraph QAOA num wires" begin
        @testset "num wires MaxKColSubgraphPhaseSeparationGate" begin
            γs = rand(length(graphs)) * 2π
            κs = [3, 2, 2]

            gates = MaxKColSubgraphPhaseSeparationGate.(γs, κs, graphs)
            for (gate, graph, κ) ∈ zip(gates, graphs, κs)
                @test Qaintessent.num_wires(gate) == graph.n * κ
            end
        end

        @testset "num wires ParityRingMixerGate" begin
            N = 3 # number of test cases
            βs = rand(N) * 2π
            ds = rand(2:9, N)

            gates = ParityRingMixerGate.(βs, ds)
            for (gate, d) ∈ zip(gates, ds)
                @test Qaintessent.num_wires(gate) == d
            end
        end

        @testset "num wires RNearbyValuesMixerGate" begin
            βs = rand(4) * 2π
            ds = [2, 3, 6, 9]
            rs = [1, 2, 1, 5]

            gates = RNearbyValuesMixerGate.(βs, rs, ds)
            for (gate, d) ∈ zip(gates, ds)
                @test Qaintessent.num_wires(gate) == d
            end
        end

        @testset "num wires PartitionMixerGate" begin
            ds = [2, 3, 6, 9]
            nums_partition_layers = [2, 3, 4, 5]

            for (d, num_partition_layers) ∈ zip(ds, nums_partition_layers)
                # create num_partition_layers partition parts, each of
                # random length and with random entries
                partition = [
                    # create h x 2 matrix with random height 2 from a random permutation of the indices...
                    (reshape(randperm(d)[1:(2*rand(1:div(d, 2)))], (:, 2))
                    # ... then map the rows to tuples
                        |> eachrow .|> Tuple)
                    for _ ∈ 1:num_partition_layers
                ]

                β = rand() * 10 + 0.1
                gate = PartitionMixerGate(β, d, partition)
                @test Qaintessent.num_wires(gate) == d
            end
        end
    end
end