using Test
using TestSetExtensions
using Qaintessent


@testset ExtendedTestSet "registers" begin
    @testset "quantum register construction" begin
        N = 3
        q1 = qreg(3)
        cgc = Circuit(q1)
        add!(X, q1(3))
        add!(Y, q1(1))
        add!(Z, q1(2))

        cgc_ref = Circuit{N}([
            circuit_gate(3, X),
            circuit_gate(1, Y),
            circuit_gate(2, Z)
        ])
        @test matrix(cgc) ≈ matrix(cgc_ref)
    end

    @testset "quantum register wire numbers with multiple registers" begin
        N = 10
        q1 = qreg(3)
        q2 = qreg(2)
        q3 = qreg(5)
        cgc = Circuit(q1,q2,q3)
        
        add!(X, q1(3))
        add!(Y, q1(1))
        add!(Z, q1(2))
        add_control!(Y, q2(1), q3(3))
        add_control!(Z, q1(2), q2(2))
        add_control!(TGate(), q3(2), q3(4))
        add_control!(SGate(), q1(2), q2(1))

        cgc_ref = Circuit{N}([
            circuit_gate(3, X),
            circuit_gate(1, Y),
            circuit_gate(2, Z),
            circuit_gate(4, Y, 8),
            circuit_gate(2, Z, 5),
            circuit_gate(7, TGate(), 9),
            circuit_gate(2, SGate(), 4),
            ])
        
        @test matrix(cgc) ≈ matrix(cgc_ref)
    end

    @testset "quantum register exceptions" begin
        N = 3
        q1 = qreg(3)

        @test_throws ErrorException("Register object has yet to be used in a Circuit object") q1[1]
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") q1(1)
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") add!(X, q1(1))
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") add!(X, q1)
    end
end