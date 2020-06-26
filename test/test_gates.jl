using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


function isunitary(g::Qaintessent.AbstractGate)
    Qaintessent.matrix(g) * Qaintessent.matrix(Base.adjoint(g)) ≈ I
end


@testset ExtendedTestSet "quantum gates" begin

    θ = 0.7*π
    ϕ = 0.4*π
    n = randn(3); n /= norm(n)

    for g in [X, Y, Z, HadamardGate(), SGate(), TGate(), RxGate(θ), RyGate(θ), RzGate(θ), RotationGate(θ, n), PhaseShiftGate(ϕ), controlled_not()]
        @test isunitary(g)
        gdag = adjoint(g)
        @test LinearAlgebra.ishermitian(g) == (Qaintessent.matrix(gdag) == Qaintessent.matrix(g))
    end

    θ = 6π
    for g in [RxGate(θ), RyGate(θ), RzGate(θ), PhaseShiftGate(θ)]
        @test isunitary(g)
        gdag = adjoint(g)
        @test LinearAlgebra.ishermitian(g) == (Qaintessent.matrix(gdag) == Qaintessent.matrix(g))
    end

    @test Qaintessent.matrix(RotationGate(θ, [1, 0, 0])) ≈ Qaintessent.matrix(RxGate(θ))
    @test Qaintessent.matrix(RotationGate(θ, [0, 1, 0])) ≈ Qaintessent.matrix(RyGate(θ))
    @test Qaintessent.matrix(RotationGate(θ, [0, 0, 1])) ≈ Qaintessent.matrix(RzGate(θ))

    #Test MatrixGate
    N = 3
    d = 2
    A = rand(ComplexF64, d^N ,d^N)
    U, R = qr(A)
    U = Array(U)
    Q = diagm(rand([-0.9:0.1:0.9...], N))
    UQ = U*Q
    GateU = MatrixGate(U)

    @test Qaintessent.matrix(GateU) ≈ U
    @test isunitary(GateU)
    @test_throws ErrorException MatrixGate(UQ)

end
