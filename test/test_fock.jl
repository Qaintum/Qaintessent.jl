using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent


# @testset ExtendedTestSet "fock state" begin
#     s = 5
#     m = 3
#     trunc = 6
#     @testset "uniform fock state constructor" begin
#         fs = FockState(s, m, trunc)
#
#         state = zeros(ComplexF64, (m, trunc))
#         state[:,s] .= 1.0 + 0.0im
#
#         @test fs.state == state
#         @test fs.pure
#     end
#
#     @testset "uniform fock state constructor exceptions" begin
#         m = -2
#         @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, m, trunc)
#
#         m = 2
#         trunc = -1
#         @test_throws ErrorException("Truncated photon number must be positive") FockState(s, m, trunc)
#
#         trunc = 2
#         @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s, m, trunc)
#     end
#
#     s = [3, 4, 1]
#     m = 3
#     trunc = 6
#     @testset "general fock state constructor" begin
#         fs = FockState(s, m, trunc)
#
#         state = zeros(ComplexF64, (m, trunc))
#         for (i,s) in enumerate(s)
#             state[i,s] = 1.0 + 0.0im
#         end
#
#         @test fs.state == state
#         @test fs.pure
#     end
#
#     @testset "uniform fock state constructor exceptions" begin
#         m = -2
#         @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, m, trunc)
#
#         m = 3
#         trunc = -1
#         @test_throws ErrorException("Truncated photon number must be positive") FockState(s, m, trunc)
#
#         trunc = 6
#         s = [1,2,3,4]
#         @test_throws ErrorException("3-mode FockState created but 4 modes provided in `s`") FockState(s, m, trunc)
#
#         trunc = 2
#         s = [3, 4, 1]
#         @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s, m, trunc)
#     end
# end

# @testset ExtendedTestSet "standard fock gates" begin
#     m = 2
#     N = 4
#     @testset "creation and annihilation operators" begin
#         creat = Qaintessent.matrix(âDag{N}())
#         dest = Qaintessent.matrix(â{N}())
#
#         @test creat * dest ≈ Qaintessent.matrix(n̂{N}())
#         @test (dest*creat)[1:N-1, 1:N-1] - (creat*dest)[1:N-1, 1:N-1] ≈ Matrix(1.0I, N-1, N-1)
#     end
#
#     @testset "phase-shift operators" begin
#         θ = rand(ComplexF64, 1)[1]
#
#         a = Qaintessent.matrix(â{N}())
#         U = Qaintessent.matrix(Û{N}(θ))
#         Udag = Qaintessent.matrix(Û{N}(-θ))
#
#         @test Udag*a*U ≈ a*exp(-im*θ)
#     end
#
# end
#
# @testset ExtendedTestSet "apply fock state" begin
#
#     s = [1,2,3]
#     m = 3
#     N = 4
#
#     fs = FockState(s, m, N)
#     @testset ExtendedTestSet "apply creation operator to fock state" begin
#         apply(âDag{N}(), fs)
#
#         @test fs.state == Complex{Float64}[0.0 1.0 0.0 0.0;
#                                            0.0 0.0  √2 0.0;
#                                            0.0 0.0 0.0  √3]
#     end
#
#     fs = FockState(s, m, N)
#     @testset ExtendedTestSet "apply destruction operator to fock state" begin
#
#         apply(â{N}(), fs)
#
#         @test fs.state == Complex{Float64}[0.0 0.0 0.0 0.0;
#                                            1.0 0.0 0.0 0.0;
#                                            0.0  √2 0.0 0.0]
#     end
#
#     fs = FockState(s, m, N)
#     @testset ExtendedTestSet "apply photon number operator to fock state" begin
#
#         apply(n̂{N}(), fs)
#
#         @test fs.state == Complex{Float64}[1.0 0.0 0.0 0.0;
#                                            0.0 2.0 0.0 0.0;
#                                            0.0 0.0 3.0 0.0]
#     end
#
#     fs = FockState(s, m, N)
#     @testset ExtendedTestSet "apply phase-shift operator to fock state" begin
#         θ = rand(Float64, 1)[1]
#         apply(Û{N}(θ), fs)
#         @test fs.state == Complex{Float64}[exp(-im*θ) 0.0 0.0 0.0;
#                                            0.0 exp(-2im*θ) 0.0 0.0;
#                                            0.0 0.0 exp(-3im*θ) 0.0]
#     end
#
#     fs = vacuum_fock_state(m, N)
#     @testset ExtendedTestSet "apply displacement operator to fock state" begin
#         α = rand(ComplexF64, 1)[1]
#         apply(D̂{N}(α), fs)
#         cs = coherent_state(α, m, N)
#         @test fs.state == cs.state
#     end
# end

@testset ExtendedTestSet "apply fock state" begin
end
