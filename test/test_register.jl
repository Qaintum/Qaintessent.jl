using Test
using TestSetExtensions
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


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
        @test sparse_matrix(cgc) ≈ sparse_matrix(cgc_ref)
    end

    @testset "quantum register array operations " begin
        N = 3
        q1 = qreg(3)

        cgc = Circuit(q1)
        
        @test all(q1[:] .== [1,2,3])
        @test all(q1[1:2] .== [1,2])
    end

    @testset "quantum register wire numbers with multiple registers" begin
        N = 12
        q1 = qreg(3)
        q2 = qreg(2)
        q3 = qreg(5)
        q4 = qreg(2)
        cgc = Circuit(q1,q2,q3,q4)
        
        add!(X, q1(3))
        add!(Y, q1(1))
        add!(Z, q1(2))
        add_control!(Y, q2(1), q3(3))
        add_control!(Z, q1(2), q2(2))
        add_control!(TGate(), q3(2), q3(4))
        add_control!(SGate(), q1(2), q2(1))
        add!(X, q1)
        add!(SwapGate(), q2)
        add!(SwapGate(), q2(1), q4(2))
        add_control!(X, q2, q4)

        cgc_ref = Circuit{N}([
            circuit_gate(3, X),
            circuit_gate(1, Y),
            circuit_gate(2, Z),
            circuit_gate(4, Y, 8),
            circuit_gate(2, Z, 5),
            circuit_gate(7, TGate(), 9),
            circuit_gate(2, SGate(), 4),
            circuit_gate(1, X),
            circuit_gate(2, X),
            circuit_gate(3, X),
            circuit_gate(4, 5, SwapGate()),
            circuit_gate(4, 12, SwapGate()),
            circuit_gate(4, X, 11),
            circuit_gate(5, X, 12)
            ])
        
        @test sparse_matrix(cgc) ≈ sparse_matrix(cgc_ref)
    end

    @testset "quantum register exceptions" begin
        N = 3
        q1 = qreg(3)
        q2 = qreg(3)
        q3 = qreg(4)
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") q1[1]
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") q1(1)
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") add!(X, q1(1))
        @test_throws ErrorException("Register object has yet to be used in a Circuit object") add!(X, q1)

        @test_throws ErrorException("Target quantum register is not included in any circuit") add_control!(X, q1, q2)
        cgc1 = Circuit(q1)
        @test_throws ErrorException("Control quantum register is not included in any circuit") add_control!(X, q1, q2)
        cgc2 = Circuit(q2)

        @test_throws ErrorException("AbstractGate SwapGate() is applied to 2 wires, quantum register of size 3 provided. Unable to determine gates to be applied.") add!(SwapGate(), q1)
        @test_throws ErrorException("Targetted Registers belong to different Circuit objects") add!(SwapGate(), q1(1), q2(1))
        @test_throws ErrorException("Target quantum register and control quantum register are used in different circuits") add_control!(X, q1, q2)
    end
end