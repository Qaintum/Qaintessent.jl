using Test
using TestSetExtensions
using Qaintessent


@testset ExtendedTestSet "test qasm" begin
    filename = "test4.qasm"
    cgc = import_file(filename)
end
