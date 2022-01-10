using Test
using TestSetExtensions
using RBNF
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------

@testset ExtendedTestSet "qasm support" begin
    @testset "qasm simple qasm2circuit" begin
        # Adding random comments and spacings in src1 to test robustness of lexer
        src = """
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
        
        @testset  "qasm simple lexer" begin
        
            tokens = Qaintessent.lex(src)
        
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
        
        @testset "qasm simple parse tokens" begin
            qasmtokens = Qaintessent.lex(src)
            rbnf_tokens = Qaintessent.parse_qasm(qasmtokens)
            rbnf_tokens_ref = Qaintessent.Struct_mainprogram(
                        RBNF.Token{:nnreal}(1, 10, 10, "2.0", 3), 
                        Any[
                            Qaintessent.Struct_decl(RBNF.Token{:id}(2, 1, 15, "qreg", 4), 
                                RBNF.Token{:id}(2, 6, 20, "a", 1), RBNF.Token{:nninteger}(2, 8, 22, "3", 1)), 
                            Qaintessent.Struct_gate(
                                Qaintessent.Struct_gatedecl(RBNF.Token{:id}(4, 6, 32, "syndrome", 8), 
                                    Qaintessent.Struct_idlist(RBNF.Token{:id}(4, 15, 41, "theta", 5), nothing), 
                                    Qaintessent.Struct_idlist(RBNF.Token{:id}(4, 22, 48, "a1", 2), 
                                    Qaintessent.Struct_idlist(RBNF.Token{:id}(4, 25, 51, "a2", 2), nothing))), 
                            Any[
                                Qaintessent.Struct_rx(RBNF.Token{:id}(6, 4, 77, "theta", 5), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(6, 11, 84, "a1", 2), nothing)), 
                                Qaintessent.Struct_iduop(RBNF.Token{:id}(7, 1, 118, "cx", 2), nothing, 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(7, 4, 121, "a1", 2), nothing), 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(7, 7, 124, "a2", 2), nothing), nothing)))]), 
                                Qaintessent.Struct_iduop(RBNF.Token{:id}(11, 1, 151, "syndrome", 8), 
                                    Qaintessent.Struct_explist(RBNF.Token{:nnreal}(11, 10, 160, "0.2", 3), nothing), 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(11, 15, 165, "a", 1), RBNF.Token{:nninteger}(11, 17, 167, "0", 1)), 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(11, 20, 170, "a", 1), RBNF.Token{:nninteger}(11, 22, 172, "2", 1)), nothing)))]
                        )
                    
            @test rbnf_tokens.ver == rbnf_tokens_ref.ver
            @test all(rbnf_tokens.prog .== rbnf_tokens_ref.prog)
        end

        @testset "qasm simple tokens to circuit" begin
            qasmtokens = Qaintessent.lex(src)
            rbnf_tokens = Qaintessent.parse_qasm(qasmtokens)
            circuit = Qaintessent.transform_qasm(rbnf_tokens)
            circuit_ref = Circuit{3}([
                    circuit_gate(1, RxGate(0.2)),
                    circuit_gate(3, X, 1)
                    ])
            @test all(circuit .≈ circuit_ref)
        end
    end


    # ##==----------------------------------------------------------------------------------------------------------------------

    @testset "qasm complex qasm2circuit" begin
        # Adding random comments and spacings in src1 to test robustness of lexer
        src = """
        OPENQASM 2.0;
        qreg a[6];
        gate syn(alpha,beta) a1,a2,a3 //gate definition
        {
        rx(alpha) a1; //rx gate with variable theta
        cx a2,a3;
        t a3;
        sdg a1;
        crz(beta) a3,a1;
        ry(beta) a2;
        }

        gate drome(gamma,delta) a1,a2,a3 //gate definition
        {
        cz a2,a3;
        s a1;
        t a3;
        cy a2,a1;
        rz(gamma) a3;
        h a2;
        crx(delta) a3,a1;
        }

        syn(sin(0.2),cos(0.5)) a[0],a[1],a[2];
        drome(exp(3),3.4*pi) a[3],a[4],a[5];
        cry(ln(5)) a[1],a[5];
        crz(exp(4)) a[2],a[4];
        h a[3];
        """

        @testset "qasm complex lexer" begin
            tokens = Qaintessent.lex(src)

            """Structure of RNBF token is:
                struct Token{T}
                    lineno :: Int64
                    colno  :: Int32
                    offset :: Int64
                    str    :: String
                    span   :: Int32
            """

            tokens_ref = RBNF.Token[
                RBNF.Token{:unnamed}(1, 1, 1, "OPENQASM", 8)
                RBNF.Token{:nnreal}(1, 10, 10, "2.0", 3)
                RBNF.Token{:unnamed}(1, 13, 13, ";", 1)
                RBNF.Token{:id}(2, 1, 15, "qreg", 4)
                RBNF.Token{:id}(2, 6, 20, "a", 1)
                RBNF.Token{:unnamed}(2, 7, 21, "[", 1)
                RBNF.Token{:nninteger}(2, 8, 22, "6", 1)
                RBNF.Token{:unnamed}(2, 9, 23, "]", 1)
                RBNF.Token{:unnamed}(2, 10, 24, ";", 1)
                RBNF.Token{:id}(3, 1, 26, "gate", 4)
                RBNF.Token{:id}(3, 6, 31, "syn", 3)
                RBNF.Token{:unnamed}(3, 9, 34, "(", 1)
                RBNF.Token{:id}(3, 10, 35, "alpha", 5)
                RBNF.Token{:unnamed}(3, 15, 40, ",", 1)
                RBNF.Token{:id}(3, 16, 41, "beta", 4)
                RBNF.Token{:unnamed}(3, 20, 45, ")", 1)
                RBNF.Token{:id}(3, 22, 47, "a1", 2)
                RBNF.Token{:unnamed}(3, 24, 49, ",", 1)
                RBNF.Token{:id}(3, 25, 50, "a2", 2)
                RBNF.Token{:unnamed}(3, 27, 52, ",", 1)
                RBNF.Token{:id}(3, 28, 53, "a3", 2)
                RBNF.Token{:unnamed}(4, 1, 74, "{", 1)
                RBNF.Token{:id}(5, 1, 76, "rx", 2)
                RBNF.Token{:unnamed}(5, 3, 78, "(", 1)
                RBNF.Token{:id}(5, 4, 79, "alpha", 5)
                RBNF.Token{:unnamed}(5, 9, 84, ")", 1)
                RBNF.Token{:id}(5, 11, 86, "a1", 2)
                RBNF.Token{:unnamed}(5, 13, 88, ";", 1)
                RBNF.Token{:id}(6, 1, 120, "cx", 2)
                RBNF.Token{:id}(6, 4, 123, "a2", 2)
                RBNF.Token{:unnamed}(6, 6, 125, ",", 1)
                RBNF.Token{:id}(6, 7, 126, "a3", 2)
                RBNF.Token{:unnamed}(6, 9, 128, ";", 1)
                RBNF.Token{:id}(7, 1, 130, "t", 1)
                RBNF.Token{:id}(7, 3, 132, "a3", 2)
                RBNF.Token{:unnamed}(7, 5, 134, ";", 1)
                RBNF.Token{:id}(8, 1, 136, "sdg", 3)
                RBNF.Token{:id}(8, 5, 140, "a1", 2)
                RBNF.Token{:unnamed}(8, 7, 142, ";", 1)
                RBNF.Token{:id}(9, 1, 144, "crz", 3)
                RBNF.Token{:unnamed}(9, 4, 147, "(", 1)
                RBNF.Token{:id}(9, 5, 148, "beta", 4)
                RBNF.Token{:unnamed}(9, 9, 152, ")", 1)
                RBNF.Token{:id}(9, 11, 154, "a3", 2)
                RBNF.Token{:unnamed}(9, 13, 156, ",", 1)
                RBNF.Token{:id}(9, 14, 157, "a1", 2)
                RBNF.Token{:unnamed}(9, 16, 159, ";", 1)
                RBNF.Token{:id}(10, 1, 161, "ry", 2)
                RBNF.Token{:unnamed}(10, 3, 163, "(", 1)
                RBNF.Token{:id}(10, 4, 164, "beta", 4)
                RBNF.Token{:unnamed}(10, 8, 168, ")", 1)
                RBNF.Token{:id}(10, 10, 170, "a2", 2)
                RBNF.Token{:unnamed}(10, 12, 172, ";", 1)
                RBNF.Token{:unnamed}(11, 1, 174, "}", 1)
                RBNF.Token{:id}(13, 1, 177, "gate", 4)
                RBNF.Token{:id}(13, 6, 182, "drome", 5)
                RBNF.Token{:unnamed}(13, 11, 187, "(", 1)
                RBNF.Token{:id}(13, 12, 188, "gamma", 5)
                RBNF.Token{:unnamed}(13, 17, 193, ",", 1)
                RBNF.Token{:id}(13, 18, 194, "delta", 5)
                RBNF.Token{:unnamed}(13, 23, 199, ")", 1)
                RBNF.Token{:id}(13, 25, 201, "a1", 2)
                RBNF.Token{:unnamed}(13, 27, 203, ",", 1)
                RBNF.Token{:id}(13, 28, 204, "a2", 2)
                RBNF.Token{:unnamed}(13, 30, 206, ",", 1)
                RBNF.Token{:id}(13, 31, 207, "a3", 2)
                RBNF.Token{:unnamed}(14, 1, 228, "{", 1)
                RBNF.Token{:id}(15, 1, 230, "cz", 2)
                RBNF.Token{:id}(15, 4, 233, "a2", 2)
                RBNF.Token{:unnamed}(15, 6, 235, ",", 1)
                RBNF.Token{:id}(15, 7, 236, "a3", 2)
                RBNF.Token{:unnamed}(15, 9, 238, ";", 1)
                RBNF.Token{:id}(16, 1, 240, "s", 1)
                RBNF.Token{:id}(16, 3, 242, "a1", 2)
                RBNF.Token{:unnamed}(16, 5, 244, ";", 1)
                RBNF.Token{:id}(17, 1, 246, "t", 1)
                RBNF.Token{:id}(17, 3, 248, "a3", 2)
                RBNF.Token{:unnamed}(17, 5, 250, ";", 1)
                RBNF.Token{:id}(18, 1, 252, "cy", 2)
                RBNF.Token{:id}(18, 4, 255, "a2", 2)
                RBNF.Token{:unnamed}(18, 6, 257, ",", 1)
                RBNF.Token{:id}(18, 7, 258, "a1", 2)
                RBNF.Token{:unnamed}(18, 9, 260, ";", 1)
                RBNF.Token{:id}(19, 1, 262, "rz", 2)
                RBNF.Token{:unnamed}(19, 3, 264, "(", 1)
                RBNF.Token{:id}(19, 4, 265, "gamma", 5)
                RBNF.Token{:unnamed}(19, 9, 270, ")", 1)
                RBNF.Token{:id}(19, 11, 272, "a3", 2)
                RBNF.Token{:unnamed}(19, 13, 274, ";", 1)
                RBNF.Token{:id}(20, 1, 276, "h", 1)
                RBNF.Token{:id}(20, 3, 278, "a2", 2)
                RBNF.Token{:unnamed}(20, 5, 280, ";", 1)
                RBNF.Token{:id}(21, 1, 282, "crx", 3)
                RBNF.Token{:unnamed}(21, 4, 285, "(", 1)
                RBNF.Token{:id}(21, 5, 286, "delta", 5)
                RBNF.Token{:unnamed}(21, 10, 291, ")", 1)
                RBNF.Token{:id}(21, 12, 293, "a3", 2)
                RBNF.Token{:unnamed}(21, 14, 295, ",", 1)
                RBNF.Token{:id}(21, 15, 296, "a1", 2)
                RBNF.Token{:unnamed}(21, 17, 298, ";", 1)
                RBNF.Token{:unnamed}(22, 1, 300, "}", 1)
                RBNF.Token{:id}(24, 1, 303, "syn", 3)
                RBNF.Token{:unnamed}(24, 4, 306, "(", 1)
                RBNF.Token{:id}(24, 5, 307, "sin", 3)
                RBNF.Token{:unnamed}(24, 8, 310, "(", 1)
                RBNF.Token{:nnreal}(24, 9, 311, "0.2", 3)
                RBNF.Token{:unnamed}(24, 12, 314, ")", 1)
                RBNF.Token{:unnamed}(24, 13, 315, ",", 1)
                RBNF.Token{:id}(24, 14, 316, "cos", 3)
                RBNF.Token{:unnamed}(24, 17, 319, "(", 1)
                RBNF.Token{:nnreal}(24, 18, 320, "0.5", 3)
                RBNF.Token{:unnamed}(24, 21, 323, ")", 1)
                RBNF.Token{:unnamed}(24, 22, 324, ")", 1)
                RBNF.Token{:id}(24, 24, 326, "a", 1)
                RBNF.Token{:unnamed}(24, 25, 327, "[", 1)
                RBNF.Token{:nninteger}(24, 26, 328, "0", 1)
                RBNF.Token{:unnamed}(24, 27, 329, "]", 1)
                RBNF.Token{:unnamed}(24, 28, 330, ",", 1)
                RBNF.Token{:id}(24, 29, 331, "a", 1)
                RBNF.Token{:unnamed}(24, 30, 332, "[", 1)
                RBNF.Token{:nninteger}(24, 31, 333, "1", 1)
                RBNF.Token{:unnamed}(24, 32, 334, "]", 1)
                RBNF.Token{:unnamed}(24, 33, 335, ",", 1)
                RBNF.Token{:id}(24, 34, 336, "a", 1)
                RBNF.Token{:unnamed}(24, 35, 337, "[", 1)
                RBNF.Token{:nninteger}(24, 36, 338, "2", 1)
                RBNF.Token{:unnamed}(24, 37, 339, "]", 1)
                RBNF.Token{:unnamed}(24, 38, 340, ";", 1)
                RBNF.Token{:id}(25, 1, 342, "drome", 5)
                RBNF.Token{:unnamed}(25, 6, 347, "(", 1)
                RBNF.Token{:id}(25, 7, 348, "exp", 3)
                RBNF.Token{:unnamed}(25, 10, 351, "(", 1)
                RBNF.Token{:nninteger}(25, 11, 352, "3", 1)
                RBNF.Token{:unnamed}(25, 12, 353, ")", 1)
                RBNF.Token{:unnamed}(25, 13, 354, ",", 1)
                RBNF.Token{:nnreal}(25, 14, 355, "3.4", 3)
                RBNF.Token{:unnamed}(25, 17, 358, "*", 1)
                RBNF.Token{:id}(25, 18, 359, "pi", 2)
                RBNF.Token{:unnamed}(25, 20, 361, ")", 1)
                RBNF.Token{:id}(25, 22, 363, "a", 1)
                RBNF.Token{:unnamed}(25, 23, 364, "[", 1)
                RBNF.Token{:nninteger}(25, 24, 365, "3", 1)
                RBNF.Token{:unnamed}(25, 25, 366, "]", 1)
                RBNF.Token{:unnamed}(25, 26, 367, ",", 1)
                RBNF.Token{:id}(25, 27, 368, "a", 1)
                RBNF.Token{:unnamed}(25, 28, 369, "[", 1)
                RBNF.Token{:nninteger}(25, 29, 370, "4", 1)
                RBNF.Token{:unnamed}(25, 30, 371, "]", 1)
                RBNF.Token{:unnamed}(25, 31, 372, ",", 1)
                RBNF.Token{:id}(25, 32, 373, "a", 1)
                RBNF.Token{:unnamed}(25, 33, 374, "[", 1)
                RBNF.Token{:nninteger}(25, 34, 375, "5", 1)
                RBNF.Token{:unnamed}(25, 35, 376, "]", 1)
                RBNF.Token{:unnamed}(25, 36, 377, ";", 1)
                RBNF.Token{:id}(26, 1, 379, "cry", 3)
                RBNF.Token{:unnamed}(26, 4, 382, "(", 1)
                RBNF.Token{:id}(26, 5, 383, "ln", 2)
                RBNF.Token{:unnamed}(26, 7, 385, "(", 1)
                RBNF.Token{:nninteger}(26, 8, 386, "5", 1)
                RBNF.Token{:unnamed}(26, 9, 387, ")", 1)
                RBNF.Token{:unnamed}(26, 10, 388, ")", 1)
                RBNF.Token{:id}(26, 12, 390, "a", 1)
                RBNF.Token{:unnamed}(26, 13, 391, "[", 1)
                RBNF.Token{:nninteger}(26, 14, 392, "1", 1)
                RBNF.Token{:unnamed}(26, 15, 393, "]", 1)
                RBNF.Token{:unnamed}(26, 16, 394, ",", 1)
                RBNF.Token{:id}(26, 17, 395, "a", 1)
                RBNF.Token{:unnamed}(26, 18, 396, "[", 1)
                RBNF.Token{:nninteger}(26, 19, 397, "5", 1)
                RBNF.Token{:unnamed}(26, 20, 398, "]", 1)
                RBNF.Token{:unnamed}(26, 21, 399, ";", 1)
                RBNF.Token{:id}(27, 1, 401, "crz", 3)
                RBNF.Token{:unnamed}(27, 4, 404, "(", 1)
                RBNF.Token{:id}(27, 5, 405, "exp", 3)
                RBNF.Token{:unnamed}(27, 8, 408, "(", 1)
                RBNF.Token{:nninteger}(27, 9, 409, "4", 1)
                RBNF.Token{:unnamed}(27, 10, 410, ")", 1)
                RBNF.Token{:unnamed}(27, 11, 411, ")", 1)
                RBNF.Token{:id}(27, 13, 413, "a", 1)
                RBNF.Token{:unnamed}(27, 14, 414, "[", 1)
                RBNF.Token{:nninteger}(27, 15, 415, "2", 1)
                RBNF.Token{:unnamed}(27, 16, 416, "]", 1)
                RBNF.Token{:unnamed}(27, 17, 417, ",", 1)
                RBNF.Token{:id}(27, 18, 418, "a", 1)
                RBNF.Token{:unnamed}(27, 19, 419, "[", 1)
                RBNF.Token{:nninteger}(27, 20, 420, "4", 1)
                RBNF.Token{:unnamed}(27, 21, 421, "]", 1)
                RBNF.Token{:unnamed}(27, 22, 422, ";", 1)
                RBNF.Token{:id}(28, 1, 424, "h", 1)
                RBNF.Token{:id}(28, 3, 426, "a", 1)
                RBNF.Token{:unnamed}(28, 4, 427, "[", 1)
                RBNF.Token{:nninteger}(28, 5, 428, "3", 1)
                RBNF.Token{:unnamed}(28, 6, 429, "]", 1)
                RBNF.Token{:unnamed}(28, 7, 430, ";", 1)
                ]
            @test all(tokens_ref .== tokens)
        end

        @testset "qasm complex parse tokens" begin
            qasmtokens = Qaintessent.lex(src)
            rbnf_tokens = Qaintessent.parse_qasm(qasmtokens)
            rbnf_tokens_ref = 
                Qaintessent.Struct_mainprogram(
                    RBNF.Token{:nnreal}(1, 10, 10, "2.0", 3), 
                    Any[
                        Qaintessent.Struct_decl(
                            RBNF.Token{:id}(2, 1, 15, "qreg", 4), 
                            RBNF.Token{:id}(2, 6, 20, "a", 1), 
                            RBNF.Token{:nninteger}(2, 8, 22, "6", 1)), 
                        Qaintessent.Struct_gate(
                            Qaintessent.Struct_gatedecl(RBNF.Token{:id}(3, 6, 31, "syn", 3), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 10, 35, "alpha", 5), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 16, 41, "beta", 4), nothing)), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 22, 47, "a1", 2), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 25, 50, "a2", 2), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(3, 28, 53, "a3", 2), nothing)))), 
                            Any[
                                Qaintessent.Struct_rx(RBNF.Token{:id}(5, 4, 79, "alpha", 5), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(5, 11, 86, "a1", 2), nothing)), 
                                Qaintessent.Struct_iduop(RBNF.Token{:id}(6, 1, 120, "cx", 2), nothing, 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(6, 4, 123, "a2", 2), nothing), 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(6, 7, 126, "a3", 2), nothing), nothing))), 
                                Qaintessent.Struct_t(Qaintessent.Struct_argument(RBNF.Token{:id}(7, 3, 132, "a3", 2), nothing)), 
                                Qaintessent.Struct_sdg(Qaintessent.Struct_argument(RBNF.Token{:id}(8, 5, 140, "a1", 2), nothing)), 
                                Qaintessent.Struct_crz(RBNF.Token{:id}(9, 5, 148, "beta", 4), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(9, 11, 154, "a3", 2), nothing), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(9, 14, 157, "a1", 2), nothing)), 
                                Qaintessent.Struct_ry(RBNF.Token{:id}(10, 4, 164, "beta", 4), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(10, 10, 170, "a2", 2), nothing))]), 
                        Qaintessent.Struct_gate(
                            Qaintessent.Struct_gatedecl(RBNF.Token{:id}(13, 6, 182, "drome", 5), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(13, 12, 188, "gamma", 5), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(13, 18, 194, "delta", 5), nothing)), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(13, 25, 201, "a1", 2), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(13, 28, 204, "a2", 2), 
                                Qaintessent.Struct_idlist(RBNF.Token{:id}(13, 31, 207, "a3", 2), nothing)))), 
                            Any[
                                Qaintessent.Struct_iduop(RBNF.Token{:id}(15, 1, 230, "cz", 2), nothing, 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(15, 4, 233, "a2", 2), nothing), 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(15, 7, 236, "a3", 2), nothing), nothing))), 
                                Qaintessent.Struct_s(Qaintessent.Struct_argument(RBNF.Token{:id}(16, 3, 242, "a1", 2), nothing)), 
                                Qaintessent.Struct_t(Qaintessent.Struct_argument(RBNF.Token{:id}(17, 3, 248, "a3", 2), nothing)), 
                                Qaintessent.Struct_iduop(RBNF.Token{:id}(18, 1, 252, "cy", 2), nothing, 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(18, 4, 255, "a2", 2), nothing), 
                                    Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(18, 7, 258, "a1", 2), nothing), nothing))), 
                                Qaintessent.Struct_rz(RBNF.Token{:id}(19, 4, 265, "gamma", 5), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(19, 11, 272, "a3", 2), nothing)), 
                                Qaintessent.Struct_h(Qaintessent.Struct_argument(RBNF.Token{:id}(20, 3, 278, "a2", 2), nothing)), 
                                Qaintessent.Struct_crx(RBNF.Token{:id}(21, 5, 286, "delta", 5), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(21, 12, 293, "a3", 2), nothing), 
                                    Qaintessent.Struct_argument(RBNF.Token{:id}(21, 15, 296, "a1", 2), nothing))]), 
                        Qaintessent.Struct_iduop(RBNF.Token{:id}(24, 1, 303, "syn", 3), 
                            Qaintessent.Struct_explist(Qaintessent.Struct_fnexp(RBNF.Token{:id}(24, 5, 307, "sin", 3), RBNF.Token{:nnreal}(24, 9, 311, "0.2", 3)), 
                            Qaintessent.Struct_explist(Qaintessent.Struct_fnexp(RBNF.Token{:id}(24, 14, 316, "cos", 3), RBNF.Token{:nnreal}(24, 18, 320, "0.5", 3)), nothing)), 
                            Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(24, 24, 326, "a", 1), RBNF.Token{:nninteger}(24, 26, 328, "0", 1)), 
                            Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(24, 29, 331, "a", 1), RBNF.Token{:nninteger}(24, 31, 333, "1", 1)), 
                            Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(24, 34, 336, "a", 1), RBNF.Token{:nninteger}(24, 36, 338, "2", 1)), nothing)))), 
                        Qaintessent.Struct_iduop(RBNF.Token{:id}(25, 1, 342, "drome", 5), 
                            Qaintessent.Struct_explist(Qaintessent.Struct_fnexp(RBNF.Token{:id}(25, 7, 348, "exp", 3), RBNF.Token{:nninteger}(25, 11, 352, "3", 1)), 
                            Qaintessent.Struct_explist(Qaintessent.Struct_bin(RBNF.Token{:nnreal}(25, 14, 355, "3.4", 3), RBNF.Token{:unnamed}(25, 17, 358, "*", 1), 
                            Qaintessent.Struct_pi()), nothing)), 
                            Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(25, 22, 363, "a", 1), RBNF.Token{:nninteger}(25, 24, 365, "3", 1)), 
                            Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(25, 27, 368, "a", 1), RBNF.Token{:nninteger}(25, 29, 370, "4", 1)), 
                            Qaintessent.Struct_mixedlist(Qaintessent.Struct_argument(RBNF.Token{:id}(25, 32, 373, "a", 1), RBNF.Token{:nninteger}(25, 34, 375, "5", 1)), nothing)))), 
                        Qaintessent.Struct_cry(Qaintessent.Struct_fnexp(RBNF.Token{:id}(26, 5, 383, "ln", 2), RBNF.Token{:nninteger}(26, 8, 386, "5", 1)), 
                            Qaintessent.Struct_argument(RBNF.Token{:id}(26, 12, 390, "a", 1), RBNF.Token{:nninteger}(26, 14, 392, "1", 1)), 
                            Qaintessent.Struct_argument(RBNF.Token{:id}(26, 17, 395, "a", 1), RBNF.Token{:nninteger}(26, 19, 397, "5", 1))), 
                        Qaintessent.Struct_crz(Qaintessent.Struct_fnexp(RBNF.Token{:id}(27, 5, 405, "exp", 3), RBNF.Token{:nninteger}(27, 9, 409, "4", 1)), 
                            Qaintessent.Struct_argument(RBNF.Token{:id}(27, 13, 413, "a", 1), RBNF.Token{:nninteger}(27, 15, 415, "2", 1)), 
                            Qaintessent.Struct_argument(RBNF.Token{:id}(27, 18, 418, "a", 1), RBNF.Token{:nninteger}(27, 20, 420, "4", 1))), 
                        Qaintessent.Struct_h(Qaintessent.Struct_argument(RBNF.Token{:id}(28, 3, 426, "a", 1), RBNF.Token{:nninteger}(28, 5, 428, "3", 1)))])

            @test rbnf_tokens.ver == rbnf_tokens_ref.ver
            @test all(rbnf_tokens.prog .== rbnf_tokens_ref.prog)
        end

        @testset "qasm complex tokens to circuit" begin
            qasmtokens = Qaintessent.lex(src)
            rbnf_tokens = Qaintessent.parse_qasm(qasmtokens)
            circuit = Qaintessent.transform_qasm(rbnf_tokens)

            α = sin(0.2)
            β = cos(0.5)
            γ = exp(3)
            δ = 3.4π
            ϵ = log(5)
            ζ = exp(4)

            circuit_ref = Circuit{6}(
                [circuit_gate(1, RxGate(α)),
                circuit_gate(3, X, 2),
                circuit_gate(3, TGate()),
                circuit_gate(1, SdagGate()),
                circuit_gate(1, RzGate(β), 3),
                circuit_gate(2, RyGate(β)),
                circuit_gate(6, Z, 5),
                circuit_gate(4, SGate()),
                circuit_gate(6, TGate()),
                circuit_gate(4, Y, 5),
                circuit_gate(6, RzGate(γ)),
                circuit_gate(5, HadamardGate()),
                circuit_gate(4, RxGate(δ), 6),
                circuit_gate(6, RyGate(ϵ), 2),
                circuit_gate(5, RzGate(ζ),3),
                circuit_gate(4, HadamardGate())]
        )
            @test all(circuit .≈ circuit_ref)
        end
    end

    # ##==----------------------------------------------------------------------------------------------------------------------


    @testset "test reading openqasm" begin
        src1 = """
        // Repetition code syndrome measurement
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg d[3];
        qreg a[2];
        qreg c[3];
        creg e[2];
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
            cx d3,a1;
            cx a1,d3;
            cx d3,a1;
        }
        syn2drome(0.1*pi) d[0],d[1],d[2],a[0],a[1];
        t d[1];
        s a[0];
        h d[0];
        x a[1];
        y d[0];
        z c[0];

        measure d[1] -> e[0];
        measure a[0] -> e[1];
        """

        d1 = qreg(3)
        a1 = qreg(2)
        c1 = qreg(3)
        cgc_ref = Circuit(d1, a1, c1)

        N = num_wires(cgc_ref)
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
            circuit_gate(4, 3, SwapGate()),
            circuit_gate(2, TGate()),
            circuit_gate(4, SGate()),
            circuit_gate(1, HadamardGate()),
            circuit_gate(5, X),
            circuit_gate(1, Y),
            circuit_gate(6, Z)
        ]
        meas_ref = [
            mop(Z, 2),
            mop(Z, 4)
        ]

        append!(cgc_ref, gates_ref)
        add_measurement!(cgc_ref, meas_ref)

        cgc = qasm2cgc(src1)    
        ψ = randn(ComplexQ, 2^N)
        
        @test apply(ψ, cgc_ref.moments) ≈ apply(ψ, cgc.moments)
        @test apply(ψ, cgc_ref) ≈ apply(ψ, cgc)
    end

    ##==----------------------------------------------------------------------------------------------------------------------

    @testset "test openqasm warnings" begin
        src1 = """
        // Repetition code syndrome measurement
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg d[3];
        qreg a[2];
        qreg c[3];
        creg e[2]
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
            cx d3,a1;
            cx a1,d3;
            cx d3,a1;
        }
        syn2drome(0.1*pi) d[0],d[1],d[2],a[0],a[1];
        t d[1];
        barrier d[1];
        s a[0];
        h d[0];

        measure d[1] -> e[0];
        measure a[0] -> e[1];
        """

        @test_logs (:warn,"Barrier operation is not yet supported") qasm2cgc(src1)
        
        src2 = """
        // Repetition code syndrome measurement
        OPENQASM 2.0;
        include "qelib1.inc";
        qreg d[3];
        qreg a[2];
        qreg c[3];
        creg e[2]
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
            cx d3,a1;
            cx a1,d3;
            cx d3,a1;
        }
        syn2drome(0.1*pi) d[0],d[1],d[2],a[0],a[1];
        t d[1];
        reset d[1];
        s a[0];
        h d[0];
        cx d[0], d[1];

        measure d[1] -> e[0];
        measure a[0] -> e[1];
        """

        @test_logs (:warn,"Reset operation is not yet supported") qasm2cgc(src2)
    end



    ##==----------------------------------------------------------------------------------------------------------------------


    @testset "qasm writeout" begin
        src_ref64 = """
        OPENQASM 2.0;

        include "qelib1.inc";

        qreg qregister[8];
        cx qregister[0],qregister[3];
        cx qregister[1],qregister[3];
        cx qregister[1],qregister[4];
        cx qregister[2],qregister[4];
        ry(0.3141592653589793) qregister[3];
        cx qregister[1],qregister[0];
        cx qregister[2],qregister[0];
        cx qregister[2],qregister[4];
        cx qregister[3],qregister[4];
        ry(0.3141592653589793) qregister[0];
        rx(0.3141592653589793) qregister[1];
        cx qregister[2],qregister[1];
        t qregister[1];
        s qregister[3];
        h qregister[0];
        tdg qregister[3];"""

        src_ref32 = """
        OPENQASM 2.0;

        include "qelib1.inc";

        qreg qregister[8];
        cx qregister[0],qregister[3];
        cx qregister[1],qregister[3];
        cx qregister[1],qregister[4];
        cx qregister[2],qregister[4];
        ry(0.31415927) qregister[3];
        cx qregister[1],qregister[0];
        cx qregister[2],qregister[0];
        cx qregister[2],qregister[4];
        cx qregister[3],qregister[4];
        ry(0.31415927) qregister[0];
        rx(0.31415927) qregister[1];
        cx qregister[2],qregister[1];
        t qregister[1];
        s qregister[3];
        h qregister[0];
        tdg qregister[3];"""


        d1 = qreg(3)
        a1 = qreg(2)
        c1 = qreg(3)
        cgc_ref = Circuit(d1, a1, c1)

        N = num_wires(cgc_ref)

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
            circuit_gate(4, TdagGate()),
        ]
        append!(cgc_ref, gates_ref)

        if FloatQ == Float32
            src_ref = src_ref32
        else
            src_ref = src_ref64
        end

        ψ = randn(ComplexQ, 2^N)
        cgc = qasm2cgc(src_ref64)

        @test apply(ψ, cgc.moments) ≈ apply(ψ, cgc_ref.moments)
        @test src_ref == cgc2qasm(cgc_ref)
    end
end