using Qaintessent

@testset ExtendedTestSet "registers" begin
    @testset "quantum register construction" begin
        N = 3
        q1 = qreg(3)
        cgc = CircuitGateChain([q1])
        cgc([
            single_qubit_circuit_gate(q1[3], X, N),
            single_qubit_circuit_gate(q1[1], Y, N),
            single_qubit_circuit_gate(q1[2], Z, N)
        ])
        cgc_ref = CircuitGateChain{N}([
            single_qubit_circuit_gate(3, X, N),
            single_qubit_circuit_gate(1, Y, N),
            single_qubit_circuit_gate(2, Z, N)
        ])
        for i in 1:3
            @test cgc[i] ≈ cgc_ref[i]
        end
    end

    @testset "quantum register wire numbers with multiple registers" begin
        N = 10
        q1 = qreg(3)
        q2 = qreg(2)
        q3 = qreg(5)
        cgc = CircuitGateChain([q1,q2,q3])
        cgc([
            single_qubit_circuit_gate(q1[3], X, N),
            single_qubit_circuit_gate(q1[1], Y, N),
            single_qubit_circuit_gate(q1[2], Z, N),
            controlled_circuit_gate(q2[1], q3[3], Y, N),
            controlled_circuit_gate(q1[2], q2[2], Z, N),
            controlled_circuit_gate(q3[2], q3[4], TGate(), N),
            controlled_circuit_gate(q1[2], q2[1], SGate(), N),

        ])
        cgc_ref = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, X, N),
        single_qubit_circuit_gate(1, Y, N),
        single_qubit_circuit_gate(2, Z, N),
        controlled_circuit_gate(4, 8, Y, N),
        controlled_circuit_gate(2, 5, Z, N),
        controlled_circuit_gate(6, 9, TGate(), N),
        controlled_circuit_gate(2, 4, SGate(), N),
        ])
        for i in 1:3
            @test cgc[i] ≈ cgc_ref[i]
        end
    end

    @testset "quantum register exceptions" begin
        N = 3
        q1 = qreg(3)

        @test_throws ErrorException("Register object has yet to be used in a CircuitGateChain object") q1[1]
    end

    @testset "classical register construction" begin
        c1 = creg(3)
        q1 = qreg(3)
        cgc = CircuitGateChain([q1], [c1])
        set_creg!(c1, 5)
        @test c1.n == [1,0,1]
        set_creg!(c1, 3)
        @test c1.n == [0,1,1]
        set_creg!(c1, 2)
        @test c1.n == [0,1,0]
    end

    @testset "classical register checks" begin
        c1 = creg(3)
        q1 = qreg(3)
        cgc = CircuitGateChain([q1], [c1])
        set_creg!(c1, 5)
        @test eval(reg_check(c1, 5))
        @test !eval(reg_check(c1, 2))
    end

    @testset "classical register exceptions" begin
        c1 = creg(3)
        q1 = qreg(2)

        @test_throws ErrorException("Unable to store integer value 10 in BitArray of size 3") set_creg!(c1, 10)

        cgc = CircuitGateChain([q1], [c1])
        @test_throws ErrorException("Unable to store integer value 14 in BitArray of size 3") set_creg!(c1, 14)
    end
end
