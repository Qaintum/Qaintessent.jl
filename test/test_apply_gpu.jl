using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent
using Qaintessent.gpu



@testset ExtendedTestSet "apply gates to statevector w/ GPU" begin
    N = 5
    ψ = rand(Qaintessent.ComplexQ, 2^N)

    @testset "apply basic gates to statevector w/ GPU" begin
        # single qubit gates
        θ = 2π*rand()
        ϕ = 2π*rand()
        n = randn(3); n /= norm(n)
        U, _ = qr(randn(ComplexF64, (2,2)))
    
        for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            i = rand(1:N)
            cg = circuit_gate(i, g)
            cga = CircuitGate{1,AbstractGate}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}
            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))
            Qaintessent.gpu.apply!(s1, cg)
            apply!(s2, cga)            
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end
    end

    
    @testset "apply multiple basic gates to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        n = randn(3); n /= norm(n)
        U, _ = qr(randn(ComplexF64, (2,2)))

        g = [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), MatrixGate(U)]
        cgs = circuit_gate.(rand(1:N, length(g)), shuffle!(g))
        cga = CircuitGate[]
        for cg in cgs
            push!(cga, CircuitGate{1,AbstractGate}(cg.iwire, cg.gate))
        end
        s1 = Statevector(deepcopy(ψ))
        s2 = Statevector(deepcopy(ψ))
        Qaintessent.gpu.apply!(s1, cgs)
        apply!(s2, cgs)
        @test s1.state ≈ s2.state
        @test s1.state ≈ sparse_matrix(cgs, N) * ψ
    end

    @testset "apply moments to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        n = randn(3); n /= norm(n)
        U, _ = qr(randn(ComplexF64, (2,2)))

        # single qubit gates
        for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ), MatrixGate(U)]
            i = rand(1:N)
            cg = Moment([circuit_gate(i, g)])
            cga = Moment([CircuitGate{1,AbstractGate}(cg[1].iwire, cg[1].gate)]) # generate same gate with type AbstractGate{1}
            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            Qaintessent.gpu.apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end
    end

    @testset "apply multiple moments to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        n = randn(3); n /= norm(n)
        U, _ = qr(randn(ComplexF64, (2,2)))

        # single qubit gates
        g = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ), MatrixGate(U)]
        i = randperm(length(g))[1:3]
        j = randperm(N)[1:3]
        cgs = circuit_gate.(j, g[i])
        cga = CircuitGate[]
        for cg in cgs
            push!(cga, CircuitGate{1,AbstractGate}(cg.iwire, cg.gate))
        end
        s1 = Statevector(deepcopy(ψ))
        s2 = Statevector(deepcopy(ψ))

        Qaintessent.gpu.apply!(s1, cgs)
        apply!(s2, Moment(cga))
        @test s1.state ≈ s2.state
        @test s1.state ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "apply basic controlled gates to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        n = randn(3); n /= norm(n)
        U, _ = qr(randn(ComplexF64, (2,2)))

        # control gate
        for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ), MatrixGate(U)]
            num_controls = rand(1:(N÷2))
            i = rand(1:N)
            j = Tuple(shuffle!([1:i-1; i+1:N][1:num_controls]))
            M = num_controls + 1
            cg = circuit_gate(i, g, j)
            cga = CircuitGate{M,ControlledGate{typeof(g)}}(cg.iwire, cg.gate)
            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            Qaintessent.gpu.apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end
    end

    @testset "apply 2-qubit gates to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        η = 2π*rand()

        for g in [SwapGate(), EntanglementXXGate(θ), EntanglementYYGate(ϕ), EntanglementZZGate(η)]
            i = rand(1:N)
            j = rand([1:i-1; i+1:N])
            cg = circuit_gate((i, j), g)
            cga = CircuitGate{2,typeof(g)}(cg.iwire, cg.gate)

            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            Qaintessent.gpu.apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end
    end

    @testset "apply multiple 2-qubit gates to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        η = 2π*rand()

        g = [SwapGate(), EntanglementXXGate(θ), EntanglementYYGate(ϕ), EntanglementZZGate(η)]
        i = randperm(length(g))[1:3]
        j = [(randperm(N)[1:2]...,) for i in 1:3]
        cgs = circuit_gate.(j, g[i])
        cga = CircuitGate[]
        for cg in cgs
            push!(cga, CircuitGate{2,AbstractGate}(cg.iwire, cg.gate))
        end

        s1 = Statevector(deepcopy(ψ))
        s2 = Statevector(deepcopy(ψ))

        Qaintessent.gpu.apply!(s1, cgs)
        apply!(s2, cga)
        @test s1.state ≈ s2.state
        @test s1.state ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "apply controlled 2-qubit gates to statevector w/ GPU" begin
        θ = 2π*rand()
        ϕ = 2π*rand()
        η = 2π*rand()
        
        for g in [SwapGate(), EntanglementXXGate(θ), EntanglementYYGate(ϕ), EntanglementZZGate(η)]
            num_controls = rand(1:(N÷2))
            i = rand(1:N)
            j = rand([1:i-1; i+1:N])
            k = Tuple(shuffle!(setdiff(1:N, [i,j]))[1:num_controls])
            M = num_controls + 2
            cg = circuit_gate((i, j), g, k)
            cga = CircuitGate{M,ControlledGate{typeof(g)}}(cg.iwire, cg.gate)

            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            Qaintessent.gpu.apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end
    end

    @testset "apply k-qubit MatrixGate to statevector w/ GPU" begin
        #  MatrixGate: k qubits
        # k = rand(1:(N÷2))
        k = 3
        A = rand(ComplexF64, 2^k, 2^k)
        U, R = qr(A)
        U = Array(U);
        g = MatrixGate(U)
        iwire = randperm(N)[1:k]
        sort!(iwire)
        cg = circuit_gate((iwire...,), g)
        cga = CircuitGate{k,AbstractGate}((iwire...,), g)
    
        m = sparse_matrix(cga, N)

        s1 = Statevector(deepcopy(ψ))
        s2 = Statevector(deepcopy(ψ))

        Qaintessent.gpu.apply!(s1, cg)
        s3 = apply(s2.state, cga)
        apply!(s2, cga)
        
        @test s1.state ≈ s2.state
        @test s1.state ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "apply controlled k-qubit MatrixGate to statevector w/ GPU" begin
        # MatrixGate: k qubits
        k = rand(1:(N÷2))
        num_controls = rand(1:(N-k))
        A = rand(ComplexF64, 2^k, 2^k)
        U, R = qr(A)
        U = Array(U);
        g = MatrixGate(U)
        iwire = Tuple(sort(randperm(N)[1:k]))
        control = Tuple(shuffle!(setdiff(1:N, iwire))[1:num_controls])
        M = num_controls + k
        cg = circuit_gate(iwire, g, control)
        cga = CircuitGate{M,ControlledGate{<:AbstractGate}}((iwire..., control...,), cg.gate)

        s1 = Statevector(deepcopy(ψ))
        s2 = Statevector(deepcopy(ψ))

        Qaintessent.gpu.apply!(s1, cg)
        s3 = apply(s2.state, cga)
        apply!(s2, cga)
        @test s1.state ≈ s2.state
        @test s1.state ≈ sparse_matrix(cga, N) * ψ
    end

    @testset "circuit measurement w/ GPU" begin
        cgs = [
            circuit_gate(rand(1:N), X),
            circuit_gate(rand(1:N), HadamardGate()),
            circuit_gate(rand(1:N), Z),
            circuit_gate(rand(1:N), Y),
            ]

        iwires = [(rand(1:N),), (rand(1:N),), (rand(1:N),)]
        meas = mop.([Z, Z, Z], iwires)
        c = Circuit{N}(cgs, meas)
        s = Statevector(deepcopy(ψ))
        s2 = Statevector(deepcopy(ψ))
        ψs = apply(ψ, cgs)
        output = Qaintessent.gpu.apply!(s, c)
        
        @test s.state ≈ ψs
        @test [dot(ψs, apply(ψs, m)) for m in meas] ≈ output
    end 

    @testset "circuit measurement with arbitrary observables w/ GPU" begin
        cgs = [
            circuit_gate(rand(1:N), X),
            circuit_gate(rand(1:N), HadamardGate()),
            circuit_gate(rand(1:N), Z),
            circuit_gate(rand(1:N), Y),
            ]

        iwires = [(rand(1:N),), (rand(1:N),), (rand(1:N),)]
        meas = mop.([matrix(Z), matrix(Z), matrix(Z)], iwires)
        c = Circuit{N}(cgs, meas)
        s = Statevector(deepcopy(ψ))
        ψs = apply(ψ, cgs)
        output = Qaintessent.gpu.apply!(s, c)

        @test s.state ≈ ψs
        @test [dot(ψs, apply(ψs, m)) for m in meas] ≈ output
    end 
end