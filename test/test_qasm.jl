using Qaintessent


@testset ExtendedTestSet "test reading openqasm" begin
    src1 = """
    // Repetition code syndrome measurement
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg d[3];
    qreg a[2];
    qreg c[3];
    creg syn[4];
    gate syndrome d1,d2,d3,a1,a2
    {
      cx d1,a1; cx d2,a1;
      cx d2,a2; cx d3,a2;
    }
    gate syn2drome(theta) d1,d2,d3,a1,a2
    {
        syndrome d1,d2,d3,a1,a2;
        syndrome d2,d3,a1,d1,a2;
        rx(theta) d2;
    }
    syn2drome(2*pi) d[0],d[1],d[2],a[0],a[1];
    if(syn==0) cx d[2],d[1];
    y d[1];
    """

    syn1 = creg(4)
    d1 = qreg(3)
    a1 = qreg(2)
    c1 = qreg(3)
    cgc_ref = CircuitGateChain([d1, a1, c1], [syn1])

    N = size(cgc_ref)

    gates_ref = [
        controlled_circuit_gate(1, 4, X, N),
        controlled_circuit_gate(2, 4, X, N),
        controlled_circuit_gate(2, 5, X, N),
        controlled_circuit_gate(3, 5, X, N),
        controlled_circuit_gate(2, 1, X, N),
        controlled_circuit_gate(3, 1, X, N),
        controlled_circuit_gate(3, 5, X, N),
        controlled_circuit_gate(4, 5, X, N),
        single_qubit_circuit_gate(2, RxGate(2π), N),
        controlled_circuit_gate(3, 2, X, N),
        single_qubit_circuit_gate(2, Y, N),
    ]
    cgc_ref(gates_ref)

    cgc = rbnfparse(src1)

    ψ = randn(ComplexF64, 2^N)
    println(exp(2))
    ans_ref = apply(cgc_ref, ψ)

    @test ans_ref ≈ ans_ref
end
