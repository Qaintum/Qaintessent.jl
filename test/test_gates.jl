using Test

using Qaintessent


@testset "quantum gates" begin

    @test Qaintessent.matrix(controlled_not()) ≈ [1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]

end
