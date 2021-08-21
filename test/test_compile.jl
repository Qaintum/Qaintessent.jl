using Test
using TestSetExtensions
using LinearAlgebra
using Random
using RandomMatrices
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------
import LinearAlgebra.BLAS: @blasfunc


using LinearAlgebra: BlasInt
for (s, elty) in (("dlarfg_", Float64),
                  ("zlarfg_", ComplexF64))
    @eval begin
        function RandomMatrices.larfg!(n::Int, α::Ptr{$elty}, x::Ptr{$elty}, incx::Int, τ::Ptr{$elty})
	    ccall((@blasfunc($s), LAPACK.liblapack), Nothing,
		  (Ptr{Int}, Ptr{$elty}, Ptr{$elty}, Ptr{Int}, Ptr{$elty}),
                  Ref(n), α, x, Ref(incx), τ)
        end
    end
end


@testset ExtendedTestSet "compilediagonal unitaries helper functions" begin
    @testset "greyencode" begin
        g = Qaintessent.greyencode.(collect(0:15))
        @test all(g .== [0, 1, 3, 2, 6, 7, 5, 4, 12, 13, 15, 14, 10, 11, 9, 8])
    end

    @testset "svalue" begin
        nref = abs(rand(Int, 1)[1])%10
        b = collect(Qaintessent.svalue(nref))
        n = 0
        for i in b
            n += 1 << (i-1)
        end
        @test n == nref
    end

    @testset "fillψ!" begin
        N = 2
        d = exp.(im .* rand(Float64, 2^N))
        l = 2^(N-1)
        ψ = zeros(Float64, l-1)
        Qaintessent.fillψ!(d, ψ, l)
        ψref = zeros(Float64, l-1)
        for i in StepRange(1,1,l-1)
            ψref[i] = imag(log(d[2i-1]*d[2i+2]/(d[2i]*d[2i+1])))
        end
        @test all(ψref .≈ ψ)
    end

    @testset "ηcol" begin
        N = 4
        ref = [[-1, 0, 0, 0], [1, -1, 0, 0], [0, 1, -1, 0], [0, 0, 1, -1], [0, 0, 0, 1]]
        for i in 0:N
            @test all(ref[i+1] .== Qaintessent.ηcol(N, i))
        end
    end

    @testset "flip_state" begin
        N = 8
        @test all([1, 0, -1, 0, 1, 0, -1] .== Qaintessent.flip_state(3, N))
    end
end


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "compile2qubit unitaries helper functions" begin
    @testset "decomposeSO4" begin
        N = 2
        E = 1/sqrt(2) .* [1 im 0 0; 0 0 im 1; 0 0 im -1; 1 -im 0 0]

        M = Stewart(Float64, 2^N)
        Q, _ = qr(M)
        Q = Matrix{ComplexF64}(Q)
        if det(Q) == -1
            Q = Q*Q
        end
        A, B = Qaintessent.decomposeSO4(Q*Q)
        @test E'*kron(A, B)*E ≈ Q*Q
    end
end


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "general compile unitaries helper functions" begin
    @testset "QR decomposition" begin
        N = 2
        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        Uref = deepcopy(U)

        @testset "unblocked QR decomposition" begin
            QR, τ = Qaintessent.qr_unblocked(deepcopy(U))
            R = diag(QR)
            Q = Matrix{ComplexF64}(I, (2^N, 2^N))
            Id = Matrix{ComplexF64}(I, (2^N, 2^N))

            for i in 1:size(U)[1]-1
                u = pushfirst!(QR[i+1:2^N, i], 1)
                u = u ./ norm(u)
                H = deepcopy(Id)
                H[i:2^N, i:2^N] = H[i:2^N, i:2^N] - 2*u*u'
                Q = H*Q
                U = H*U
            end
            d = diag(U)

            @test diagm(d) ≈ U
            @test diag(U) ≈ R
            @test Q*Uref ≈ U
            @test inv(Q)*U ≈ Uref
        end

        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        Uref = deepcopy(U)

        @testset "unblocked! QR decomposition" begin
            QR = deepcopy(U)
            τ = Qaintessent.qr_unblocked!(QR)
            R = diag(QR)
            Q = Matrix{ComplexF64}(I, (2^N, 2^N))
            Id = Matrix{ComplexF64}(I, (2^N, 2^N))

            for i in 1:size(U)[1]-1
                u = pushfirst!(QR[i+1:2^N, i], 1)
                u = u ./ norm(u)
                H = deepcopy(Id)
                H[i:2^N, i:2^N] = H[i:2^N, i:2^N] - 2*u*u'
                Q = H*Q
                U = H*U
            end
            d = diag(U)

            @test diagm(d) ≈ U
            @test diag(U) ≈ R
            @test Q*Uref ≈ U
            @test inv(Q)*U ≈ Uref
        end
    end


    @testset "stateprep" begin
        @inline function allequal(x)
            length(x) < 2 && return true
            e1 = x[1]
            i = 2
            @inbounds for i=2:length(x)
                x[i] == e1 || return false
            end
            return true
        end

        N = 3
        ψ = rand(ComplexF64, 2^N)
        ψ = ψ ./ norm(ψ)
        angles = exp.(im.*angle.(ψ))
        ψ = ψ ./ angles
        ψref = deepcopy(ψ)
        ϕ = zeros(ComplexF64, (2^N))
        ϕ[1] = 1

        for j in 1:N-1
            cg = Qaintessent.stateprep(ψ[1:2^(j-1):end], N, j)
            ψ = apply(ψ, cg)
        end
        θ = real(atan(-ψ[2^(N-1)+1]./ψ[1]).*2)

        if !isnan(θ)
            cg = CircuitGate((N,), RyGate(θ))
            ψ = apply(ψ, cg)
        end

        @test ψ ≈ ϕ
    end

    @testset "inverseM" begin
        N = 4
        mref = [ 0.125   0.125   0.125   0.125   0.125   0.125   0.125   0.125;
                0.125  -0.125   0.125  -0.125   0.125  -0.125   0.125  -0.125;
                0.125  -0.125  -0.125   0.125   0.125  -0.125  -0.125   0.125;
                0.125   0.125  -0.125  -0.125   0.125   0.125  -0.125  -0.125;
                0.125   0.125  -0.125  -0.125  -0.125  -0.125   0.125   0.125;
                0.125  -0.125  -0.125   0.125  -0.125   0.125   0.125  -0.125;
                0.125  -0.125   0.125  -0.125  -0.125   0.125  -0.125   0.125;
                0.125   0.125   0.125   0.125  -0.125  -0.125  -0.125  -0.125
                ]
        @test Qaintessent.inverseM(N) ≈ mref
    end
end


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "unitary compilation" begin
    @testset "general compile" begin
        N = 6
        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        M = Stewart(ComplexF64, 2^N)

        cgs = unitary2circuit(deepcopy(U), N)

        ψ = rand(ComplexF64, 2^N)
        ψ_ref = U*ψ
        ψ_compiled = apply(ψ, cgs)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "general compile exceptions" begin
        N = 6
        U = rand(ComplexF64, (2^N, 2^N))
        @test_throws ErrorException("Only unitary matrices can be compiled into a quantum circuit") unitary2circuit(U)

        U = rand(ComplexF64, (2^N, 2^(N-1)))
        @test_throws ErrorException("Only unitary matrices can be compiled into a quantum circuit") unitary2circuit(U)
    end

    @testset "compile 1 qubit standard unitary" begin
        N = 1
        M = Stewart(ComplexF64, 2)
        random_θ = rand(Float64, 3)
        for gate in [X, Y, Z, HadamardGate(), TGate(), SGate(), RxGate(random_θ[1]), RyGate(random_θ[2]), RzGate(random_θ[3])]
            U = matrix(gate)

            cgs = unitary2circuit(U)

            ψ = rand(ComplexF64, 2^N)

            ψ_ref = U*ψ
            ψ_compiled = apply(ψ, cgs)
            @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
        end
    end

    @testset "compile 1 qubit random unitary" begin
        N = 1
        U, _ = qr(Stewart(ComplexF64, 2^N))
        U = Matrix(U)
        M = Stewart(ComplexF64, 2)

        cgs = unitary2circuit(U, N)

        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(ψ, cgs)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "compile diagonal unitaries" begin
        N = 5
        U = diagm(exp.(im .* rand(Float64, 2^N)))
        M = Stewart(ComplexF64, 2^N)

        cgs = unitary2circuit(U, N)

        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(ψ, cgs)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "compile 2 qubit unitaries" begin
        N = 2
        U, _ = qr(Matrix(Stewart(ComplexF64, 2^N)))
        U = Matrix(U)
        M = Stewart(ComplexF64, 2^N)

        cgs = unitary2circuit(U, N)
        ψ = rand(ComplexF64, 2^N)

        ψ_ref = U*ψ
        ψ_compiled = apply(ψ, cgs)

        @test ψ_ref'*M*ψ_ref ≈ ψ_compiled'*M*ψ_compiled
    end

    @testset "compile standard 2 qubit unitaries" begin
        N = 2
        M = Stewart(ComplexF64, 2^N)
        l = deepcopy(M)
        for _ in 1:400
            random_θ = rand(Float64, 3)
            for gate in [X, Y, Z, TGate(), SGate(), RxGate(random_θ[1]), RyGate(random_θ[2]), RzGate(random_θ[3])]
                
                nU = Matrix(sparse_matrix(circuit_gate(1, gate, 2)))
                cgs = unitary2circuit(nU, N)

                ψ = rand(ComplexF64, 2^N)
                ψ_ref = nU*ψ
                ψ_compiled = apply(ψ, cgs)
                
                @test isapprox(ψ_ref'*(M*ψ_ref), ψ_compiled'*(M*ψ_compiled), rtol=1e-5, atol=1e-5)
            end
        end
    end
end
