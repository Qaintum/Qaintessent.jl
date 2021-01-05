using Test
using TestSetExtensions
using RBNF
using Qaintessent


@testset ExtendedTestSet "test lexer" begin
    # Adding random comments and spacings in src1 to test robustness of lexer
    src1 = """
    OPENQASM 2.0;
    qreg a[3];

    gate syndrome(theta) a1,a2 //gate definition
    {
    rx(theta) a1; //rx gate with variable theta
    cx a1,a2;
    //controlled x gate
    }

    syndrome(0.2) a[0],a[2];
    """

    tokens = Qaintessent.lex(src1)

    """Structure of RNBF token is:
        struct Token{T}
            lineno :: Int64
            colno  :: Int32
            offset :: Int64
            str    :: String
            span   :: Int32
    """

    tokens_ref = RBNF.Token[
        # OPENQASM starting lines
        RBNF.Token{:unnamed}(1, 1, 1, "OPENQASM", 8),
        RBNF.Token{:nnreal}(1, 10, 10, "2.0", 3),
        RBNF.Token{:unnamed}(1, 13, 13, ";", 1),
        # Setting up quantum registers
        RBNF.Token{:id}(2, 1, 15, "qreg", 4),
        RBNF.Token{:id}(2, 6, 20, "a", 1),
        RBNF.Token{:unnamed}(2, 7, 21, "[", 1),
        RBNF.Token{:nninteger}(2, 8, 22, "3", 1),
        RBNF.Token{:unnamed}(2, 9, 23, "]", 1),
        RBNF.Token{:unnamed}(2, 10, 24, ";", 1),
        # Seting up custom gate
        RBNF.Token{:id}(4, 1, 27, "gate", 4),
        RBNF.Token{:id}(4, 6, 32, "syndrome", 8),
        RBNF.Token{:unnamed}(4, 14, 40, "(", 1),
        RBNF.Token{:id}(4, 15, 41, "theta", 5),
        RBNF.Token{:unnamed}(4, 20, 46, ")", 1),
        RBNF.Token{:id}(4, 22, 48, "a1", 2),
        RBNF.Token{:unnamed}(4, 24, 50, ",", 1),
        RBNF.Token{:id}(4, 25, 51, "a2", 2),
        RBNF.Token{:unnamed}(5, 1, 72, "{", 1),
        RBNF.Token{:id}(6, 1, 74, "rx", 2),
        RBNF.Token{:unnamed}(6, 3, 76, "(", 1),
        RBNF.Token{:id}(6, 4, 77, "theta", 5),
        RBNF.Token{:unnamed}(6, 9, 82, ")", 1),
        RBNF.Token{:id}(6, 11, 84, "a1", 2),
        RBNF.Token{:unnamed}(6, 13, 86, ";", 1),
        RBNF.Token{:id}(7, 1, 118, "cx", 2),
        RBNF.Token{:id}(7, 4, 121, "a1", 2),
        RBNF.Token{:unnamed}(7, 6, 123, ",", 1),
        RBNF.Token{:id}(7, 7, 124, "a2", 2),
        RBNF.Token{:unnamed}(7, 9, 126, ";", 1),
        RBNF.Token{:unnamed}(9, 1, 148, "}", 1),
        # Using customg ate
        RBNF.Token{:id}(11, 1, 151, "syndrome", 8),
        RBNF.Token{:unnamed}(11, 9, 159, "(", 1),
        RBNF.Token{:nnreal}(11, 10, 160, "0.2", 3),
        RBNF.Token{:unnamed}(11, 13, 163, ")", 1),
        RBNF.Token{:id}(11, 15, 165, "a", 1),
        RBNF.Token{:unnamed}(11, 16, 166, "[", 1),
        RBNF.Token{:nninteger}(11, 17, 167, "0", 1),
        RBNF.Token{:unnamed}(11, 18, 168, "]", 1),
        RBNF.Token{:unnamed}(11, 19, 169, ",", 1),
        RBNF.Token{:id}(11, 20, 170, "a", 1),
        RBNF.Token{:unnamed}(11, 21, 171, "[", 1),
        RBNF.Token{:nninteger}(11, 22, 172, "2", 1),
        RBNF.Token{:unnamed}(11, 23, 173, "]", 1),
        RBNF.Token{:unnamed}(11, 24, 174, ";", 1)
    ]
    @test all(tokens_ref .== tokens)
end

@testset ExtendedTestSet "test parse_qasm" begin
    src2 = """
    OPENQASM 2.0;
    qreg a[3];
    gate syndrome(theta) a1,a2
    {
    rx(theta) a1;
    cx a1,a2;
    }
    syndrome(0.2) a[0],a[2];
    """

    qasmtokens = Qaintessent.lex(src2)
    rbnf_tokens = Qaintessent.parse_qasm(qasmtokens)
    rbnf_tokens_ref = Qaintessent.Struct_mainprogram(
            RBNF.Token{:nnreal}(1, 10, 10, "2.0", 3),
            Any[
                Qaintessent.Struct_decl(
                    RBNF.Token{:id}(2, 1, 15, "qreg", 4),
                    RBNF.Token{:id}(2, 6, 20, "a", 1),
                    RBNF.Token{:nninteger}(2, 8, 22, "3", 1)),
                Qaintessent.Struct_gate(
                    Qaintessent.Struct_gatedecl(
                        RBNF.Token{:id}(3, 6, 31, "syndrome", 8),
                        Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 15, 40, "theta", 5), nothing),
                        Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 22, 47, "a1", 2),
                        Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 25, 50, "a2", 2), nothing))),
                        Any[Qaintessent.Struct_rx(RBNF.Token{:id}(5, 4, 58, "theta", 5),
                        Qaintessent.Struct_argument(RBNF.Token{:id}(5, 11, 65, "a1", 2), nothing)),
                        Qaintessent.Struct_iduop(RBNF.Token{:id}(6, 1, 69, "cx", 2), nothing,
                        Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(6, 4, 72, "a1", 2), nothing),
                        Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(6, 7, 75, "a2", 2), nothing), nothing)))]),
                Qaintessent.Struct_iduop(RBNF.Token{:id}(8, 1, 81, "syndrome", 8),
                Qaintessent.Struct_explist(RBNF.Token{:nnreal}(8, 10, 90, "0.2", 3), nothing),
                Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(8, 15, 95, "a", 1),
                    RBNF.Token{:nninteger}(8, 17, 97, "0", 1)),
                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(8, 20, 100, "a", 1),
                    RBNF.Token{:nninteger}(8, 22, 102, "2", 1)), nothing)))
                ]
            )
        # @test rbnf_tokens_ref.program .== rbnf_tokens.program
        @test rbnf_tokens == rbnf_tokens
end


@testset ExtendedTestSet "test reading openqasm" begin
    src1 = """
    // Repetition code syndrome measurement
    OPENQASM 2.0;
    include "qelib1.inc";
    qreg d[3];
    qreg a[2];
    qreg c[3];
    gate syndrome(alpha) d1,d2,d3,a1,a2
    {
      cx d1,a1; cx d2,a1;
      cx d2,a2; cx d3,a2; ry(alpha) a1;
    }
    gate syn2drome(theta) d1,d2,d3,a1,a2
    {
        syndrome(theta) d1,d2,d3,a1,a2;
        syndrome(theta) d2,d3,a1,d1,a2;
        rx(theta) d2;
    }
    syn2drome(0.1*pi) d[0],d[1],d[2],a[0],a[1];
    t d[1];
    s a[0];
    h d[0];
    """

    d1 = qreg(3)
    a1 = qreg(2)
    c1 = qreg(3)
    cgc_ref = Circuit(d1, a1, c1)

    N = size(cgc_ref)

    gates_ref = [
        circuit_gate(4, X, 1),
        circuit_gate(4, X, 2),
        circuit_gate(5, X, 2),
        circuit_gate(5, X, 3),
        circuit_gate(4, RyGate(0.1π)),
        circuit_gate(1, X, 2),
        circuit_gate(1, X, 3),
        circuit_gate(5, X, 3),
        circuit_gate(5, X, 4),
        circuit_gate(1, RyGate(0.1π)),
        circuit_gate(2, RxGate(0.1π)),        
        circuit_gate(2, TGate()),
        circuit_gate(4, SGate()),
        circuit_gate(1, HadamardGate()),
    ]

    append!(cgc_ref, gates_ref)

    cgc = qasm2cgc(src1)

    ψ = randn(ComplexF64, 2^N)
    
    @test apply(cgc_ref.moments, ψ) ≈ apply(cgc.moments, ψ)
end


@testset ExtendedTestSet "test writing" begin
    src_ref = """
    OPENQASM 2.0;

    include "qelib1.inc";

    qreg qregister[8];
    cx qregister[1],qregister[4];
    cx qregister[2],qregister[4];
    cx qregister[2],qregister[5];
    cx qregister[3],qregister[5];
    ry(0.3141592653589793) qregister[4];
    cx qregister[2],qregister[1];
    cx qregister[3],qregister[1];
    cx qregister[3],qregister[5];
    cx qregister[4],qregister[5];
    ry(0.3141592653589793) qregister[1];
    rx(0.3141592653589793) qregister[2];
    cx qregister[3],qregister[2];
    t qregister[2];
    s qregister[4];
    h qregister[1];"""

    d1 = qreg(3)
    a1 = qreg(2)
    c1 = qreg(3)
    cgc_ref = Circuit(d1, a1, c1)

    N = size(cgc_ref)

    gates_ref = [
        circuit_gate(4, X, 1),
        circuit_gate(4, X, 2),
        circuit_gate(5, X, 2),
        circuit_gate(5, X, 3),
        circuit_gate(4, RyGate(0.1π)),
        circuit_gate(1, X, 2),
        circuit_gate(1, X, 3),
        circuit_gate(5, X, 3),
        circuit_gate(5, X, 4),
        circuit_gate(1, RyGate(0.1π)),
        circuit_gate(2, RxGate(0.1π)),
        circuit_gate(2, X, 3),
        circuit_gate(2, TGate()),
        circuit_gate(4, SGate()),
        circuit_gate(1, HadamardGate()),
    ]
    append!(cgc_ref, gates_ref)

    @test src_ref == cgc2qasm(cgc_ref)
end
