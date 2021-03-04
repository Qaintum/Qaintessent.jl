using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------

@testset ExtendedTestSet "measurement operator construction" begin
    N = 5
    meas1 = MeasurementOperator.([3 * Matrix{Float64}(I, 2^N, 2^N), 2 * Matrix{Float64}(I, 2^N, 2^N)], ((1,2,3,4,5),))
    @test Qaintessent.check_commute(meas1) == true

    meas2 = MeasurementOperator.([X, Y], [(1,),(2,)])
    @test Qaintessent.check_commute(meas2) == true
end

@testset ExtendedTestSet "measurement operator abstract matrix exceptions" begin
    N = 2
    small_matrix = Matrix{Float64}(I, 2^(N)-1, 2^(N)-1)
    large_matrix = Matrix{Float64}(I, 2^(N)+1, 2^(N)+1)
    uneven_matrix = Matrix{Float64}(I, 2^(N), 2^(N+1))
    iwires = (3,4)
    non_hermitian_matrix = Qaintessent.sparse_matrix(circuit_gate(1, RxGate(0.2), 2))

    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(small_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(small_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(large_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(large_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(uneven_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be a $(2^N) × $(2^N) matrix.") MeasurementOperator(uneven_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be Hermitian.") MeasurementOperator(non_hermitian_matrix, iwires)
    @test_throws ErrorException("Measurement operator must be Hermitian.") MeasurementOperator(non_hermitian_matrix, iwires)

    noncomm_matrix_1 = sparse_matrix(circuit_gate(1, HadamardGate(), 2))
    noncomm_matrix_2 = sparse_matrix(circuit_gate(1, X, 2))
    m = MeasurementOperator.([noncomm_matrix_1, noncomm_matrix_2], (iwires,))
    @test Qaintessent.check_commute(m) == false
end

@testset ExtendedTestSet "measurement operator abstract matrix exceptions" begin
    iwires = (1,)

    non_hermitian_matrix = RxGate(0.2)
    @test_throws ErrorException("Measurement operator must be Hermitian.") MeasurementOperator(non_hermitian_matrix, iwires)
    iwires = (1,)
    noncomm_gate_1 = HadamardGate()
    noncomm_gate_2 = X
    m = MeasurementOperator.([noncomm_gate_1, noncomm_gate_2], (iwires,))
    @test Qaintessent.check_commute(m) == false
end

@testset ExtendedTestSet "circuit construction" begin
    N = 1
    cgs = [
        circuit_gate(1, X),
        circuit_gate(1, HadamardGate()),
        circuit_gate(1, Z),
        circuit_gate(1, Y),
        ]
    ref_moments = Moment.(cgs)
    meas = MeasurementOperator(Matrix{Float64}(I, 2^N, 2^N), (1,))
    c = Circuit{N}(cgs, [meas])
    for i in 1:length(cgs)
        @test c[i] ≈ ref_moments[i]
    end
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.moments, ψ)

    ψs = apply(cgs, ψ)
    ψref = apply(c, ψ)
    
    @test [dot(ψs, apply(meas, ψs))] ≈ apply(c, ψ)
end

@testset ExtendedTestSet "circuit construction with circuit gate measurements" begin
    N = 3
    cgc = [
        circuit_gate(1, X),
        circuit_gate(2, HadamardGate()),
        circuit_gate(3, Z),
        circuit_gate(1, Y),
    ]
    iwires = [(1,), (2,), (3,)]
    meas = MeasurementOperator.([X, X, X], iwires)
    c = Circuit{N}(cgc, meas)
    ψ = randn(ComplexF64, 2^N)
    @test distribution(c, ψ) ≈ apply(c.moments, ψ)

    ψs = apply(cgc, ψ)
    @test [dot(ψs, apply(m, ψs)) for m in meas] ≈ apply(c, ψ)
end


@testset ExtendedTestSet "reduced density matrix" begin
    N = 4
    ψ = randn(ComplexF64, 2^N)
    χ = randn(ComplexF64, 2^N)
    # full density matrix |ψ⟩⟨χ|
    ρ = reshape(kron(conj(χ), ψ), 2^N, 2^N)

    id = Matrix{ComplexF64}(I, 2, 2)
    A = randn(ComplexF64, 2, 2)
    B = randn(ComplexF64, 2, 2)
    @testset ExtendedTestSet "reduced density matrix correctness" begin
        @test sum(kron(A, B) .* rdm(N, (4, 2), ψ, χ)) ≈ sum(kron(B, id, A, id) .* ρ)
    end

    @testset ExtendedTestSet "reduced density matrix exceptions" begin
        @test_throws ErrorException("Need at least one wire to act on.") rdm(N, (), ψ, χ)
        @test_throws ErrorException("Number of gate wires cannot be larger than total number of wires.") rdm(N, (1, 2, 3, 4, 5), ψ, χ)
        @test_throws ErrorException("Wire indices must be unique.") rdm(N, (2, 2), ψ, χ)
        @test_throws ErrorException("Wire index cannot be smaller than 1.") rdm(N, (-1, 2), ψ, χ)
        @test_throws ErrorException("Wire index cannot be larger than total number of wires.") rdm(N, (5, 1), ψ, χ)
    end
end
