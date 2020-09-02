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
