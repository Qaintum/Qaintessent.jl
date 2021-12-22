using Test
using TestSetExtensions
using LinearAlgebra
using Random
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------
@testset ExtendedTestSet "apply gates to statevector" begin
    @testset ExtendedTestSet "apply gates to statevector helper functions" begin
        @testset "test flipqubit!" begin
            N = 5
            n = rand(1:N)
            ψ = Statevector(N)
            ref = zeros(Int, 2^N)
            for i in 1:2^N
                bits = digits(i-1, base=2, pad=N)
                bits[n] = 1 - bits[n]
                ref[i] = sum(bits .* 2 .^[0:N-1...]) + 1
            end

            Qaintessent.flipqubit!(ψ, n)
            @test all(ψ.perm .== ref)
        end

        @testset "test 1 qubit orderqubit!" begin
            N = 5
            n = rand(1:N)
            ψ = Statevector(N)
            ref = zeros(Int, 2^N)
            for i in 1:2^N
                bits = digits(i-1, base=2, pad=N)
                temp = bits[n]
                bits[n:end-1] .= bits[n+1:end]
                bits[end] = temp
                ref[sum(bits .* 2 .^[0:N-1...]) + 1] = i
            end
            Qaintessent.orderqubit!(ψ, (n,))
            @test all(ψ.perm .== ref)
        end

        @testset "test 1 qubit orderqubit! with control" begin
            N = 5
            c,n = randperm(N)[1:2]
            other_wires = setdiff(1:N, [c,n])
            ψ = Statevector(N)
            ref = zeros(Int, 2^N)
            for i in 1:2^N
                bits = digits(i-1, base=2, pad=N)
                tempc = bits[c]
                temp1 = bits[n]
                bits[1:end-2] = bits[other_wires]
                bits[end-1] = temp1
                bits[end] = tempc
                ref[sum(bits .* 2 .^[0:N-1...]) + 1] = i
            end
            Qaintessent.orderqubit!(ψ, (n,), (c,))
            @test all(ψ.perm .== ref)
        end


        @testset "test 2 qubit orderqubit!" begin
            N = 5
            n1, n2 = sort(randperm(N)[1:2])
            other_wires = setdiff(1:N, [n1,n2])
            ψ = Statevector(N)
            
            ref = zeros(Int, 2^N)
            for i in 1:2^N
                bits = digits(i-1, base=2, pad=N)
                temp1 = bits[n1]
                temp2 = bits[n2]
                bits[1:end-2] = bits[other_wires]
                bits[end-1] = temp1
                bits[end] = temp2
                ref[sum(bits .* 2 .^[0:N-1...]) + 1] = i
            end

            Qaintessent.orderqubit!(ψ, (n1, n2))
            @test all(ψ.perm .== ref)
        end

        @testset "test 2 qubit orderqubit! with control" begin
            N = 5
            wires = randperm(N)
            c = wires[1]
            n1, n2 = sort(wires[2:3])
            other_wires = setdiff(1:N, [c, n1, n2])
            ψ = Statevector(N)
            
            ref = zeros(Int, 2^N)
            for i in 1:2^N
                bits = digits(i-1, base=2, pad=N)
                temp1 = bits[n1]
                temp2 = bits[n2]
                tempc = bits[c]
                bits[1:end-3] = bits[other_wires]
                bits[end-2] = temp1
                bits[end-1] = temp2
                bits[end] = tempc
                ref[sum(bits .* 2 .^[0:N-1...]) + 1] = i
            end

            Qaintessent.orderqubit!(ψ, (n1, n2), (c,))
            @test all(ψ.perm .== ref)
        end


        @testset "test swapqubit!" begin
            N = 5
            n1 = rand(1:N)
            n2 = rand(setdiff(1:N, n1))
            ψ = Statevector(N)

            ref = zeros(Int, 2^N)
            for i in 1:2^N
                bits = digits(i-1, base=2, pad=N)
                temp = bits[n1]
                bits[n1] = bits[n2]
                bits[n2] = temp
                ref[i] = sum(bits .* 2 .^[0:N-1...]) + 1
            end

            Qaintessent.swapqubit!(ψ, n1, n2)
            @test all(ψ.perm .== ref)
        end
    end

    @testset ExtendedTestSet "apply gates to statevectors" begin
        N = 3
        θ = 0.7π
        ϕ = 0.4π
        n = randn(3); n /= norm(n)
        ψ = rand(Qaintessent.ComplexQ, 2^N)
        
        @testset "apply basic gates to statevector" begin
            # single qubit gates
            for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
                i = rand(1:N)
                cg = circuit_gate(i, g)
                cga = CircuitGate{1,AbstractGate}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}
                s1 = Statevector(deepcopy(ψ))
                s2 = Statevector(deepcopy(ψ))
                apply!(s1, cg)
                apply!(s2, cga)            
                @test s1.state ≈ s2.state
                @test s1.state ≈ sparse_matrix(cga, N) * ψ
            end
        end

        @testset "apply multiple basic gates to statevector" begin
            # single qubit gates
            g = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            i = randperm(N)[1:3]
            j = randperm(N)[1:3]
            cgs = circuit_gate.(j, shuffle!(g[i]))
            cga = CircuitGate[]
            for cg in cgs
                push!(cga, CircuitGate{1,AbstractGate}(cg.iwire, cg.gate))
            end

            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))
            apply!(s1, cgs)
            apply!(s2, cga)            
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end

        @testset "apply moments to statevector" begin
            # single qubit gates
            for g in [X, Y, Z, HadamardGate(), SGate(), SdagGate(), TGate(), TdagGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
                i = rand(1:N)
                cg = Moment([circuit_gate(i, g)])
                cga = Moment([CircuitGate{1,AbstractGate}(cg[1].iwire, cg[1].gate)]) # generate same gate with type AbstractGate{1}
                s1 = Statevector(deepcopy(ψ))
                s2 = Statevector(deepcopy(ψ))

                apply!(s1, cg)
                apply!(s2, cga)
                @test s1.state ≈ s2.state
                @test s1.state ≈ sparse_matrix(cga, N) * ψ
            end
        end

        @testset "apply multiple moments to statevector" begin
            # single qubit gates
            g = [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
            i = randperm(N)[1:3]
            j = randperm(N)[1:3]
            cgs = circuit_gate.(j, g[i])
            cga = CircuitGate[]
            for cg in cgs
                push!(cga, CircuitGate{1,AbstractGate}(cg.iwire, cg.gate))
            end
            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            apply!(s1, Moment(cgs))
            apply!(s2, Moment(cga))
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end

        @testset "apply basic controlled gates to statevector" begin
            # control gate
            for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ)]
                i = rand(1:N)
                j = rand([1:i-1; i+1:N])
                cg = circuit_gate(i, g, j)
                cga = CircuitGate{2,AbstractGate}(cg.iwire, cg.gate)
                s1 = Statevector(deepcopy(ψ))
                s2 = Statevector(deepcopy(ψ))

                apply!(s1, cg)
                apply!(s2, cga)
                @test sparse_matrix(cga, N) ≈ sparse_matrix(cg, N)
                @test s1.state ≈ s2.state
                @test s1.state ≈ sparse_matrix(cga, N) * ψ
            end
        end

        @testset "apply swap gate to statevector" begin
            θ = 2π*rand()
            ϕ = 2π*rand()
            η = 2π*rand()

            for g in [SwapGate(), EntanglementXXGate(θ), EntanglementYYGate(ϕ), EntanglementZZGate(η)]
                i = rand(1:N)
                j = rand([1:i-1; i+1:N])
                cg = circuit_gate((i, j), SwapGate())
                cga = CircuitGate{2,SwapGate}(cg.iwire, cg.gate)

                s1 = Statevector(deepcopy(ψ))
                s2 = Statevector(deepcopy(ψ))

                apply!(s1, cg)
                apply!(s2, cga)
                @test s1.state ≈ s2.state
                @test s1.state ≈ sparse_matrix(cga, N) * ψ
            end
        end

        @testset "apply controlled swap gate to statevector" begin
            θ = 2π*rand()
            ϕ = 2π*rand()
            η = 2π*rand()
            for g in [SwapGate(), EntanglementXXGate(θ), EntanglementYYGate(ϕ), EntanglementZZGate(η)]
                i = rand(1:N)
                j = rand([1:i-1; i+1:N])
                k = rand(setdiff(1:N, [i,j]))
                cg = circuit_gate((i, j), SwapGate(), k)
                cga = CircuitGate{3,ControlledGate{SwapGate}}(cg.iwire, cg.gate)

                s1 = Statevector(deepcopy(ψ))
                s2 = Statevector(deepcopy(ψ))

                apply!(s1, cg)
                apply!(s2, cga)
                @test s1.state ≈ s2.state
                @test s1.state ≈ sparse_matrix(cga, N) * ψ
            end
        end

        @testset "apply 1-qubit MatrixGate to statevector" begin
            # MatrixGate: one qubit
            d = 2
            A = rand(ComplexQ, d, d)
            U, R = qr(A)
            U = Array(U);
            g = MatrixGate(U)
            i = rand(1:N)
            cg = circuit_gate(i, g)
            cga = CircuitGate{1,MatrixGate}(cg.iwire, cg.gate) # generate same gate with type AbstractGate{1}

            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ

        end

        @testset "apply k-qubit MatrixGate to statevector" begin
            # MatrixGate: k qubits
            k = rand(1:N)
            A = rand(ComplexQ, 2^k, 2^k)
            U, R = qr(A)
            U = Array(U);
            g = MatrixGate(U)
            iwire = [rand(1:N)]
            for j in 1:k-1
                l = rand(1:N-j)
                i = setdiff([1:N...], iwire)[l]
                push!(iwire, i)
            end
            sort!(iwire)
            cg = circuit_gate((iwire...,), g)
            cga = CircuitGate{k,AbstractGate}((iwire...,), g)
            m = sparse_matrix(cga, N)

            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end

        @testset "apply controlled k-qubit MatrixGate to statevector" begin
            # MatrixGate: k qubits
            k = rand(1:N-2)
            A = rand(ComplexQ, 2^k, 2^k)
            U, R = qr(A)
            U = Array(U);
            g = MatrixGate(U)
            iwire = [rand(1:N)]
            for j in 1:k-1
                l = rand(1:N-j)
                i = setdiff([1:N...], iwire)[l]
                push!(iwire, i)
            end

            cntrl = (rand(setdiff([1:N...], iwire)),)
            
            sort!(iwire)
            cg = circuit_gate((iwire...,), g, cntrl)
            cga = CircuitGate{k+1,ControlledGate{typeof(g)}}(cg.iwire, cg.gate)
            m = sparse_matrix(cga, N)

            s1 = Statevector(deepcopy(ψ))
            s2 = Statevector(deepcopy(ψ))

            apply!(s1, cg)
            apply!(s2, cga)
            @test s1.state ≈ s2.state
            @test s1.state ≈ sparse_matrix(cga, N) * ψ
        end
    end
end
