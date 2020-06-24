using Test
using TestSetExtensions


@testset "All the tests" begin
    @includetests ["test_qasm"]
    # @includetests ARGS
end
