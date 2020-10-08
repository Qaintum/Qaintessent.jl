using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
using PyPlot


@testset ExtendedTestSet "fock state" begin
    s = 3
    m = 4
    N = 5

    @testset "uniform fock state constructor" begin
        fs = FockState(s, m, N)

        fs_ref = FockState{m,N}(ComplexF64[
                                    0 0 0 1.0 0 0;
                                    0 0 0 1.0 0 0;
                                    0 0 0 1.0 0 0;
                                    0 0 0 1.0 0 0], true)
        @test fs ≈ fs_ref
    end

    @testset "uniform fock state constructor exceptions" begin
        neg_m = -2
        @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, neg_m, N)

        neg_N = -1
        @test_throws ErrorException("Truncated photon number must be positive") FockState(s, m, neg_N)

        small_N = 2
        @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s, m, small_N)
    end

    s = [1, 2, 3, 4]
    @testset "fock state photon number as array" begin

        fs = FockState(s, m, N)

        fs_ref = FockState{m,N}(ComplexF64[
                                    0 1.0 0 0 0 0;
                                    0 0 1.0 0 0 0;
                                    0 0 0 1.0 0 0;
                                    0 0 0 0 1.0 0], true)
        @test fs ≈ fs_ref
    end

    @testset "uniform fock state constructor exceptions" begin
        neg_m = -2
        @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, neg_m, N)

        neg_N = -1
        @test_throws ErrorException("Truncated photon number must be positive") FockState(s, m, neg_N)

        diff_N = 6
        s_wrong_m = [1, 2, 3, 4, 2]
        @test_throws ErrorException("4-mode FockState created but 5 modes provided in `s`") FockState(s_wrong_m, m, diff_N)

        s_large = [3, 8, 1, 2]
        @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_large, m, N)
    end
end

@testset ExtendedTestSet "standard fock gates" begin
    m = 3
    N = 4
    @testset "creation and annihilation operators" begin
        creat = Qaintessent.matrix(âDag{N}())
        dest = Qaintessent.matrix(â{N}())

        @test creat * dest ≈ Qaintessent.matrix(n̂{N}())
        @test (dest*creat)[1:N-1, 1:N-1] - (creat*dest)[1:N-1, 1:N-1] ≈ Matrix(1.0I, N-1, N-1)


        s = [1,2,3]
        fs = FockState(s, m, N)

        cg = â{N}()
        fs2 = apply(cg, deepcopy(fs); normalize=false)
        fs_ref = FockState{m,N}(ComplexF64[
                                    1 0 0 0 0;
                                    0 sqrt(2) 0 0 0;
                                    0 0 sqrt(3) 0 0], true)
        @test fs2 ≈ fs_ref

        cg = âDag{N}()
        fs2 = apply(cg, deepcopy(fs); normalize=false)
        fs_ref = FockState{m,N}(ComplexF64[
                                    0 0 sqrt(2) 0 0;
                                    0 0 0 sqrt(3) 0;
                                    0 0 0 0 sqrt(4)], true)
        @test fs2 ≈ fs_ref
    end

    @testset "photon number operator" begin
        s = [4,1,3]
        fs = FockState(s, m, N)

        cg = n̂{N}()
        fs2 = apply(cg, deepcopy(fs); normalize=false)
        fs_ref = FockState{m,N}(ComplexF64[
                                    0 0 0 0 4;
                                    0 1 0 0 0;
                                    0 0 0 3 0], true)
        @test fs2 ≈ fs_ref
    end

    @testset "phase-shift operators" begin
        θ = rand(ComplexF64, 1)[1]

        a = Qaintessent.matrix(â{N}())
        U = Qaintessent.matrix(Û{N}(θ))
        Udag = Qaintessent.matrix(Û{N}(-θ))

        @test Udag*a*U ≈ a*exp(-im*θ)
    end

    @testset "squeezing operators" begin
        m = 1
        N = 6
        sq = Ŝ{N}(2)
        # fs = FockState(s, m, N)
        fs = vacuum_fock_state(m, N)
        apply(sq, fs; normalize=false)

        xvec = vec([-5:0.1:5...])
        yvec = vec([-5:0.1:5...])
        w = wigner(fs, xvec, yvec)
        contour(xvec, yvec, w)
    end
end

@testset ExtendedTestSet "apply fock state" begin

    s = [1, 2, 3]
    m = 3
    N = 4


    @testset "apply creation operator to fock state" begin
        fs = FockState(s, m, N)
        apply(âDag{N}(), fs; normalize=false)
        @test fs.state == Complex{Float64}[0.0 0.0  √2 0.0 0.0;
                                           0.0 0.0 0.0  √3 0.0;
                                           0.0 0.0 0.0 0.0   2]
    end


    @testset "apply destruction operator to fock state" begin
        fs = FockState(s, m, N)
        apply(â{N}(), fs; normalize=false)

        @test fs.state ≈ Complex{Float64}[1.0 0.0 0.0 0.0 0.0;
                                          0.0  √2 0.0 0.0 0.0;
                                          0.0 0.0  √3 0.0 0.0]
    end

    @testset "apply photon number operator to fock state" begin
        fs = FockState(s, m, N)
        apply(n̂{N}(), fs; normalize=false)

        @test fs.state ≈ Complex{Float64}[0.0 1.0 0.0 0.0 0.0;
                                          0.0 0.0 2.0 0.0 0.0;
                                          0.0 0.0 0.0 3.0 0.0]
    end

    @testset "apply phase-shift operator to fock state" begin
        fs = FockState(s, m, N)
        θ = rand(Float64, 1)[1]
        apply(Û{N}(θ), fs; normalize=false)

        @test fs.state ≈ Complex{Float64}[0.0 exp(-im*θ) 0.0 0.0 0.0;
                                          0.0 0.0 exp(-2im*θ) 0.0 0.0;
                                          0.0 0.0 0.0 exp(-3im*θ) 0.0]
    end

    N = 16
    @testset "apply displacement operator to fock state" begin
        fs = vacuum_fock_state(m, N)
        α = rand(ComplexF64, 1)[1]
        α /= norm(α)

        apply(D̂{N}(α), fs; normalize=false)

        cs = coherent_state(α, m, N)

        @test fs.state ≈ cs.state
    end
end

@testset ExtendedTestSet "coherent states" begin
    m = 1
    N = 15
    @testset "coherent states" begin
        α = rand(ComplexF64, 1)[1]
        cs = coherent_state(α, m, N)
        ref_norm = norm(cs.state)
        apply(â{N}(), cs; normalize=false)
        @test norm(cs.state./α) ≈ ref_norm
    end

    @testset "phase-shift operators on coherent states" begin
        α = rand(ComplexF64, 1)[1]
        θ = rand(ComplexF64, 1)[1]

        cs = coherent_state(α, m, N)
        apply(Û{N}(θ), cs; normalize=false)

        cs_shifted = coherent_state(α*exp(-im*θ), m, N)
        @test cs.state/norm(cs.state) ≈ cs_shifted.state/norm(cs_shifted.state)
    end
end

@testset ExtendedTestSet "fock state photon number as integer" begin
    modes = 4
    trunc = 5
    s = 3
    fs = FockState(s, modes, trunc)

    fs_ref = FockState{modes,trunc}(ComplexF64[
                                0 0 1.0 0 0;
                                0 0 1.0 0 0;
                                0 0 1.0 0 0;
                                0 0 1.0 0 0], true)
    @test fs ≈ fs_ref
end

@testset ExtendedTestSet "fock state photon number as integer exceptions" begin
    modes = 4
    trunc = 5
    s = 2

    neg_modes = -1
    @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, neg_modes, trunc)

    neg_trunc = -1
    @test_throws ErrorException("Truncated photon number must be positive") FockState(s, modes, neg_trunc)

    s_neg = -5
    @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)

    s_large = 7
    @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)
end

@testset ExtendedTestSet "fock state photon number as array" begin
    modes = 4
    trunc = 5
    s = [1,2,3,4]
    fs = FockState(s, modes, trunc)

    fs_ref = FockState{modes,trunc}(ComplexF64[
                                1.0 0 0 0 0;
                                0 1.0 0 0 0;
                                0 0 1.0 0 0;
                                0 0 0 1.0 0], true)
    @test fs ≈ fs_ref
end

@testset ExtendedTestSet "fock state exceptions" begin
    modes = 4
    trunc = 5
    s = [1,2,3,4]

    neg_modes = -1
    @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, neg_modes, trunc)

    neg_trunc = -1
    @test_throws ErrorException("Truncated photon number must be positive") FockState(s, modes, neg_trunc)

    s_short = [1,2,4]
    @test_throws ErrorException("4-mode FockState created but 3 modes provided in `s`") FockState(s_short, modes, trunc)

    s_neg = [-1,2,3,4]
    @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)

    s_large = [8,2,3,4]
    @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)
end

@testset ExtendedTestSet "annihilation and creation operators" begin
    modes = 3
    trunc = 4
    s = [1,2,3]
    fs = FockState(s, modes, trunc)

    cg = â{trunc}()
    fs2 = apply(cg, deepcopy(fs))
    fs_ref = FockState{modes,trunc}(ComplexF64[
                                0 0 0 0;
                                1.0 0 0 0;
                                0 sqrt(2) 0 0], true)
    @test fs2 ≈ fs_ref

    cg = âDag{trunc}()
    fs2 = apply(cg, deepcopy(fs))
    fs_ref = FockState{modes,trunc}(ComplexF64[
                                0 1.0 0 0;
                                0 0 sqrt(2) 0;
                                0 0 0 sqrt(3)], true)
    @test fs2 ≈ fs_ref
end

@testset ExtendedTestSet "photon number operator" begin
    modes = 3
    trunc = 4
    s = [4,1,3]
    fs = FockState(s, modes, trunc)

    cg = n̂{trunc}()
    fs2 = apply(cg, deepcopy(fs))
    fs_ref = FockState{modes,trunc}(ComplexF64[
                                0 0 0 4;
                                1 0 0 0;
                                0 0 3 0], true)
    @test fs2 ≈ fs_ref
end

@testset ExtendedTestSet "coherent states" begin
    modes = 3
    trunc = 8

    # alpha = (0.3+0.2im)/norm(0.3+0.2im)
    alpha = randn(ComplexF64, 1)[1]
    # alpha = alpha/norm(alpha)

    cg1 = D̂{trunc}(alpha)
    cg3 = D̂{trunc}(-alpha)
    cg2 = â{trunc}()

    # println(norm(D1))
    D = Qaintessent.matrix(cg1)
    Ddag = Qaintessent.matrix(cg3)
    # println(norm(D2))
    a = Qaintessent.matrix(cg2)

    ψ = zeros(ComplexF64, (trunc,))
    ψ[1] = 1
    ψ = D*ψ
    println([norm(a)^2 for a in ψ])
    println(sum([norm(a)^2 for a in ψ]))
    ψ2 = a*ψ

    α = ψ2./ψ
    β = α[1]
    println(β)
    α = α/β
    println([norm(a) for a in α])
end
#
@testset ExtendedTestSet "phase-shift operator" begin
    modes = 3
    trunc = 5

    theta = 0.34

    cg1 = D̂{trunc}(theta)
    cg2 = â{trunc}()

    D1 = adjoint(Qaintessent.matrix(cg1))
    # println(norm(D1))
    D2 = Qaintessent.matrix(cg1)
    # println(norm(D2))
    a = Qaintessent.matrix(cg2)

    ψ = zeros(ComplexF64, (trunc,))
    ψ[1] = 1
    ψ = D2*ψ
    ψ2 = [abs(a^2) for a in ψ]
    println(ψ2)
    println(norm(D2*ψ))
end

@testset ExtendedTestSet "displacement operator" begin
    modes = 3
    trunc = 5

    alpha = 1im

    cg1 = D̂{trunc}(alpha)
    cg2 = â{trunc}()

    D1 = adjoint(Qaintessent.matrix(cg1))
    # println(norm(D1))
    D2 = Qaintessent.matrix(cg1)
    # println(norm(D2))
    a = Qaintessent.matrix(cg2)

    ψ = [1, 0, 0, 0, 0]
    ψ = D2*ψ
    ψ2 = [abs(a^2) for a in ψ]
    println(ψ2)
    println(norm(D2*ψ))
end

@testset ExtendedTestSet "single mode squeezing operator" begin
    modes = 3
    trunc = 5

    alpha = 1im

    cg1 = D̂{trunc}(alpha)
    cg2 = â{trunc}()

    D1 = adjoint(Qaintessent.matrix(cg1))
    # println(norm(D1))
    D2 = Qaintessent.matrix(cg1)
    # println(norm(D2))
    a = Qaintessent.matrix(cg2)

    ψ = [1, 0, 0, 0, 0]
    ψ = D2*ψ
    ψ2 = [abs(a^2) for a in ψ]
    println(ψ2)
    println(norm(D2*ψ))
end

@testset ExtendedTestSet "two mode squeezing operator" begin
    modes = 3
    trunc = 5

    alpha = 1im

    cg1 = D̂{trunc}(alpha)
    cg2 = â{trunc}()

    D1 = adjoint(Qaintessent.matrix(cg1))
    # println(norm(D1))
    D2 = Qaintessent.matrix(cg1)
    # println(norm(D2))
    a = Qaintessent.matrix(cg2)

    ψ = [1, 0, 0, 0, 0]
    ψ = D2*ψ
    ψ2 = [abs(a^2) for a in ψ]
    println(ψ2)
    println(norm(D2*ψ))
end
