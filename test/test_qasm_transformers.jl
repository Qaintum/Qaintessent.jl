using Test
using TestSetExtensions
using RBNF
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------

@testset ExtendedTestSet "test lexer" begin
    # Adding random comments and spacings in src1 to test robustness of lexer
    src1 = """
    OPENQASM 2.0;
    qreg a[3];

    gate syndrome(theta) a1,a2 //gate definition
    {
    rx(theta) a1; //rx gate with variable theta
    cx a1,a2;
    t a2;
    tdg a1;
    //controlled x gate
    }

    syndrome(0.2) a[0],a[2];
    x a[2];
    cy a[1], a[2];
    z a[0];
    h a[0];
    cz a[1], a[2];
    s a[1];
    y a[2];
    rz(sin(0.3)) a[0];
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
        RBNF.Token{:unnamed}(1, 1, 1, "OPENQASM", 8)
        RBNF.Token{:nnreal}(1, 10, 10, "2.0", 3)
        RBNF.Token{:unnamed}(1, 13, 13, ";", 1)
        RBNF.Token{:id}(2, 1, 15, "qreg", 4)
        RBNF.Token{:id}(2, 6, 20, "a", 1)
        RBNF.Token{:unnamed}(2, 7, 21, "[", 1)
        RBNF.Token{:nninteger}(2, 8, 22, "3", 1)
        RBNF.Token{:unnamed}(2, 9, 23, "]", 1)
        RBNF.Token{:unnamed}(2, 10, 24, ";", 1)
        RBNF.Token{:id}(4, 1, 27, "gate", 4)
        RBNF.Token{:id}(4, 6, 32, "syndrome", 8)
        RBNF.Token{:unnamed}(4, 14, 40, "(", 1)
        RBNF.Token{:id}(4, 15, 41, "theta", 5)
        RBNF.Token{:unnamed}(4, 20, 46, ")", 1)
        RBNF.Token{:id}(4, 22, 48, "a1", 2)
        RBNF.Token{:unnamed}(4, 24, 50, ",", 1)
        RBNF.Token{:id}(4, 25, 51, "a2", 2)
        RBNF.Token{:unnamed}(5, 1, 72, "{", 1)
        RBNF.Token{:id}(6, 1, 74, "rx", 2)
        RBNF.Token{:unnamed}(6, 3, 76, "(", 1)
        RBNF.Token{:id}(6, 4, 77, "theta", 5)
        RBNF.Token{:unnamed}(6, 9, 82, ")", 1)
        RBNF.Token{:id}(6, 11, 84, "a1", 2)
        RBNF.Token{:unnamed}(6, 13, 86, ";", 1)
        RBNF.Token{:id}(7, 1, 118, "cx", 2)
        RBNF.Token{:id}(7, 4, 121, "a1", 2)
        RBNF.Token{:unnamed}(7, 6, 123, ",", 1)
        RBNF.Token{:id}(7, 7, 124, "a2", 2)
        RBNF.Token{:unnamed}(7, 9, 126, ";", 1)
        RBNF.Token{:id}(8, 1, 128, "t", 1)
        RBNF.Token{:id}(8, 3, 130, "a2", 2)
        RBNF.Token{:unnamed}(8, 5, 132, ";", 1)
        RBNF.Token{:id}(9, 1, 134, "tdg", 3)
        RBNF.Token{:id}(9, 5, 138, "a1", 2)
        RBNF.Token{:unnamed}(9, 7, 140, ";", 1)
        RBNF.Token{:unnamed}(11, 1, 162, "}", 1)
        RBNF.Token{:id}(13, 1, 165, "syndrome", 8)
        RBNF.Token{:unnamed}(13, 9, 173, "(", 1)
        RBNF.Token{:nnreal}(13, 10, 174, "0.2", 3)
        RBNF.Token{:unnamed}(13, 13, 177, ")", 1)
        RBNF.Token{:id}(13, 15, 179, "a", 1)
        RBNF.Token{:unnamed}(13, 16, 180, "[", 1)
        RBNF.Token{:nninteger}(13, 17, 181, "0", 1)
        RBNF.Token{:unnamed}(13, 18, 182, "]", 1)
        RBNF.Token{:unnamed}(13, 19, 183, ",", 1)
        RBNF.Token{:id}(13, 20, 184, "a", 1)
        RBNF.Token{:unnamed}(13, 21, 185, "[", 1)
        RBNF.Token{:nninteger}(13, 22, 186, "2", 1)
        RBNF.Token{:unnamed}(13, 23, 187, "]", 1)
        RBNF.Token{:unnamed}(13, 24, 188, ";", 1)
        RBNF.Token{:id}(14, 1, 190, "x", 1)
        RBNF.Token{:id}(14, 3, 192, "a", 1)
        RBNF.Token{:unnamed}(14, 4, 193, "[", 1)
        RBNF.Token{:nninteger}(14, 5, 194, "2", 1)
        RBNF.Token{:unnamed}(14, 6, 195, "]", 1)
        RBNF.Token{:unnamed}(14, 7, 196, ";", 1)
        RBNF.Token{:id}(15, 1, 198, "cy", 2)
        RBNF.Token{:id}(15, 4, 201, "a", 1)
        RBNF.Token{:unnamed}(15, 5, 202, "[", 1)
        RBNF.Token{:nninteger}(15, 6, 203, "1", 1)
        RBNF.Token{:unnamed}(15, 7, 204, "]", 1)
        RBNF.Token{:unnamed}(15, 8, 205, ",", 1)
        RBNF.Token{:id}(15, 10, 207, "a", 1)
        RBNF.Token{:unnamed}(15, 11, 208, "[", 1)
        RBNF.Token{:nninteger}(15, 12, 209, "2", 1)
        RBNF.Token{:unnamed}(15, 13, 210, "]", 1)
        RBNF.Token{:unnamed}(15, 14, 211, ";", 1)
        RBNF.Token{:id}(16, 1, 213, "z", 1)
        RBNF.Token{:id}(16, 3, 215, "a", 1)
        RBNF.Token{:unnamed}(16, 4, 216, "[", 1)
        RBNF.Token{:nninteger}(16, 5, 217, "0", 1)
        RBNF.Token{:unnamed}(16, 6, 218, "]", 1)
        RBNF.Token{:unnamed}(16, 7, 219, ";", 1)
        RBNF.Token{:id}(17, 1, 221, "h", 1)
        RBNF.Token{:id}(17, 3, 223, "a", 1)
        RBNF.Token{:unnamed}(17, 4, 224, "[", 1)
        RBNF.Token{:nninteger}(17, 5, 225, "0", 1)
        RBNF.Token{:unnamed}(17, 6, 226, "]", 1)
        RBNF.Token{:unnamed}(17, 7, 227, ";", 1)
        RBNF.Token{:id}(18, 1, 229, "cz", 2)
        RBNF.Token{:id}(18, 4, 232, "a", 1)
        RBNF.Token{:unnamed}(18, 5, 233, "[", 1)
        RBNF.Token{:nninteger}(18, 6, 234, "1", 1)
        RBNF.Token{:unnamed}(18, 7, 235, "]", 1)
        RBNF.Token{:unnamed}(18, 8, 236, ",", 1)
        RBNF.Token{:id}(18, 10, 238, "a", 1)
        RBNF.Token{:unnamed}(18, 11, 239, "[", 1)
        RBNF.Token{:nninteger}(18, 12, 240, "2", 1)
        RBNF.Token{:unnamed}(18, 13, 241, "]", 1)
        RBNF.Token{:unnamed}(18, 14, 242, ";", 1)
        RBNF.Token{:id}(19, 1, 244, "s", 1)
        RBNF.Token{:id}(19, 3, 246, "a", 1)
        RBNF.Token{:unnamed}(19, 4, 247, "[", 1)
        RBNF.Token{:nninteger}(19, 5, 248, "1", 1)
        RBNF.Token{:unnamed}(19, 6, 249, "]", 1)
        RBNF.Token{:unnamed}(19, 7, 250, ";", 1)
        RBNF.Token{:id}(20, 1, 252, "y", 1)
        RBNF.Token{:id}(20, 3, 254, "a", 1)
        RBNF.Token{:unnamed}(20, 4, 255, "[", 1)
        RBNF.Token{:nninteger}(20, 5, 256, "2", 1)
        RBNF.Token{:unnamed}(20, 6, 257, "]", 1)
        RBNF.Token{:unnamed}(20, 7, 258, ";", 1)
        RBNF.Token{:id}(21, 1, 260, "rz", 2)
        RBNF.Token{:unnamed}(21, 3, 262, "(", 1)
        RBNF.Token{:id}(21, 4, 263, "sin", 3)
        RBNF.Token{:unnamed}(21, 7, 266, "(", 1)
        RBNF.Token{:nnreal}(21, 8, 267, "0.3", 3)
        RBNF.Token{:unnamed}(21, 11, 270, ")", 1)
        RBNF.Token{:unnamed}(21, 12, 271, ")", 1)
        RBNF.Token{:id}(21, 14, 273, "a", 1)
        RBNF.Token{:unnamed}(21, 15, 274, "[", 1)
        RBNF.Token{:nninteger}(21, 16, 275, "0", 1)
        RBNF.Token{:unnamed}(21, 17, 276, "]", 1)
        RBNF.Token{:unnamed}(21, 18, 277, ";", 1)
            ]
    @test all(tokens_ref .== tokens)
end
