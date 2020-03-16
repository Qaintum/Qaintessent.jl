using Test
using TestSetExtensions
using Qaintessent


@testset ExtendedTestSet "quantum gates" begin
    I = [1 0; 0 1]
    @test Qaintessent.matrix(controlled_not()) ≈ [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

    @testset "Pauli X gate" begin
        @test Qaintessent.matrix(XGate()) ≈ [0.  1.; 1.  0.]
        @test Qaintessent.matrix(XGate())*Qaintessent.matrix(Base.adjoint(XGate())) ≈ I
    end

    @testset "Pauli Y gate" begin
        @test Qaintessent.matrix(YGate()) ≈ [0. -im; im  0.]
        @test Qaintessent.matrix(YGate())*Qaintessent.matrix(Base.adjoint(YGate())) ≈ I
    end

    @testset "Pauli Z gate" begin
        @test Qaintessent.matrix(ZGate()) ≈ [1.  0.; 0. -1.]
        @test Qaintessent.matrix(ZGate())*Qaintessent.matrix(Base.adjoint(ZGate())) ≈ I
    end

    @testset "Hadamard" begin
        @test Qaintessent.matrix(HadamardGate()) ≈ [1 1; 1 -1] / sqrt(2)
        @test Qaintessent.matrix(HadamardGate())*Qaintessent.matrix(Base.adjoint(HadamardGate())) ≈ I
    end

    @testset "SGate" begin
        @test Qaintessent.matrix(SGate())*Qaintessent.matrix(Base.adjoint(SGate())) ≈ I
        @test Qaintessent.matrix(SGate())*Qaintessent.matrix(SdagGate()) ≈ I
    end

    @testset "TGate" begin
        @test Qaintessent.matrix(TGate())*Qaintessent.matrix(Base.adjoint(TGate())) ≈ I
        @test Qaintessent.matrix(TGate())*Qaintessent.matrix(TdagGate()) ≈ I
    end

    @testset "Rx gate" begin
        θ = 0.5*π
        c = cos(θ/2)
        s = sin(θ/2)
        @test Qaintessent.matrix(RxGate(θ)) ≈ [c -im*s; -im*s c]
        @test Qaintessent.matrix(RxGate(θ))*Qaintessent.matrix(Base.adjoint(RxGate(θ))) ≈ I
    end

    @testset "Rx gate" begin
        θ = 0.7*π
        c = cos(θ/2)
        s = sin(θ/2)
        @test Qaintessent.matrix(RxGate(θ)) ≈ [c -im*s; -im*s c]
        @test Qaintessent.matrix(RxGate(θ))*Qaintessent.matrix(Base.adjoint(RxGate(θ))) ≈ I
    end

    @testset "Ry gate" begin
        θ = 0.7*π
        c = cos(θ/2)
        s = sin(θ/2)
        @test Qaintessent.matrix(RyGate(θ)) ≈ [c -s; s c]
        @test Qaintessent.matrix(RyGate(θ))*Qaintessent.matrix(Base.adjoint(RyGate(θ))) ≈ I
    end

    @testset "Rz gate" begin
        θ = 0.7*π
        @test Qaintessent.matrix(RzGate(θ)) ≈ [exp(-im*θ/2) 0; 0 exp(im*θ/2)]
        @test Qaintessent.matrix(RzGate(θ))*Qaintessent.matrix(Base.adjoint(RzGate(θ))) ≈ I
    end

    @testset "Rϕ gate" begin
        ϕ = 0.7*π
        @test Qaintessent.matrix(RϕGate(ϕ)) ≈ [1 0; 0 exp(im*ϕ)]
        @test Qaintessent.matrix(RϕGate(ϕ))*Qaintessent.matrix(RϕGate(-ϕ)) ≈ I
    end

end
