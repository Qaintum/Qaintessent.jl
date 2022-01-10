using Test
using TestSetExtensions
using CUDA
using Qaintessent


@testset "All the tests" begin
    skip_tests = ["test_apply_gpu.jl", "runtests.jl"]
    if false
        @includetests ARGS
    else
        term = r"(.*\.jl)"
        matches = match.((term,), readdir())
        matches = getindex.(matches[matches.!=nothing], (1,))
        filter!(s -> s ∉ skip_tests, matches)
        for test in matches
            print(splitext(test)[1], ": ")
            Base.include(Qaintessent, test)
            println()
        end
    end
end
