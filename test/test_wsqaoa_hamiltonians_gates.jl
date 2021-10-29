using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using Random
using Qaintessent.QAOAHelperDataStructs
using Qaintessent.MaxCutWSQAOA

graphs = [
    Graph(3, [(1,2), (2,3)]),
    Graph(5, [(1,2), (2,3), (3,4), (4,5), (5,1), (1,3)]),
    Graph(4, [(1,2), (1,3), (1,4), (2,3), (2,4), (3,4)])
]

function cut_size(graph::Graph, partitioning::Vector{Int})
    length(partitioning) == graph.n || throw(ArgumentError("Length of partitioning must equal the number of vertices."))
    all(partition -> partition ⊆ 0:1, partitioning) || throw(ArgumentError("Some vertices are in an invalid partition"))
    return count(partitioning[a] != partitioning[b] for (a, b) ∈ graph.edges)
end

# Create a state ψ that corresponds to a given partitioning
function ψ_from_partitioning(n::Int, partitioning::Vector{Int})::Vector{ComplexF64}
    (n > 0) || throw(DomainError("Parameter n must be positive integers."))
    partitioning ⊆ 0:1 || throw(ArgumentError("Parameter `partitioning` may only contain values in the range 0:1."))

    # Create ψ with |1> entries corresponding to the given partitioning
    ψ = kron((1 == partitioning[vertex] ? [0.0im, 1] : [1, 0.0im] for vertex ∈ 1:n)...)
    ψ
end

@testset ExtendedTestSet "maxcut QAOA Hamiltonians and gate matrices" begin
    @testset ExtendedTestSet "maxcut subgraph QAOA Hamiltonians" begin
        # Test that the phase separation Hamiltonian used in MaxCutPhaseSeparationGate
        # implements the objective function correctly.
        @testset "Hamiltonian MaxCutPhaseSeparationGate" begin
            hamiltonians = max_cut_phase_separation_hamiltonian.(graphs)
            partitionings = [rand(0:1, graph.n) for graph ∈ graphs]
            for (H, graph, partitioning) ∈ zip(hamiltonians, graphs, partitionings)
                c_size = cut_size(graph, partitioning)
                # ψ_partitioning should be scaled for comparison
                ψ_partitioning = ψ_from_partitioning(graph.n, partitioning)
                @test H * ψ_partitioning ≈ c_size * ψ_partitioning 
            end
        end

        # Test that WS-QAOA Mixer correctly maps QAOA MaxCut 
        # Mixer(X Operator) at regularization 0.5 with the initial state |0>
        @testset "Hamiltonian WSQAOAMixerGate" begin
            ε = 0.5
            init_state_randomized = false
            for graph ∈ graphs
                partitioning = zeros(graph.n)
                gate = WSQAOAMixerGate(0., partitioning, ε, init_state_randomized)
                xMixerGate = RxMixerGate(0., graph.n)
                H = wsqaoa_mixer_hamiltonian(gate)
                RX = rx_mixer_hamiltonian(xMixerGate)
                @test H ≈ RX
            end
        end
    end

    @testset ExtendedTestSet "maxcut gates adjoints" begin
        # Test that `MaxCutPhaseSeparationGate` has the correct adjoint.
        @testset "adjoint MaxKColSubgraphPhaseSeparationGate" begin
            γs = rand(length(graphs)) * 2π

            gates = MaxCutPhaseSeparationGate.(γs, graphs)
            for i ∈ 1:length(gates)
                @test Qaintessent.matrix(gates[i]) * Qaintessent.matrix(Base.adjoint(gates[i])) ≈ I
            end
        end

        # Test that WSQAOAMixerGate has the correct adjoint
        @testset "adjoint WSQAOAMixerGate" begin
            γ = 0.0
            ε = 0.0
            init_state_randomized = false

            for i ∈ 1:length(graphs)
                partitioning = zeros(Float64, graphs[i].n)
                gate = WSQAOAMixerGate(γ, partitioning, ε, init_state_randomized)
                @test Qaintessent.matrix(gate) * Qaintessent.matrix(Base.adjoint(gate)) ≈ I
            end
        end
    end

    @testset ExtendedTestSet "maxcut WS-QAOA num wires" begin
        # Test that `MaxCutPhaseSeparationGate` has the correct number of wires
        @testset "num wires MaxCutPhaseSeparationGate" begin
            γ = 0.0
            gates = MaxCutPhaseSeparationGate.(γ, graphs)
            for (gate, graph) ∈ zip(gates, graphs)
                @test Qaintessent.num_wires(gate) == graph.n
            end
        end

        # Test that `WSQAOAMixerGate` has the correct number of wires
        @testset "num wires WSQAOAMixerGate" begin
            γ = 0.0
            partitionings = zeros.(Float64, graphs[i].n for i ∈ 1:length(graphs))
            ε = 0.0
            init_state_randomized = false

            gates = WSQAOAMixerGate.(γ, partitionings, ε, init_state_randomized)
            for (gate, graph) ∈ zip(gates, graphs)
                @test Qaintessent.num_wires(gate) == graph.n
            end
        end
    end
end
