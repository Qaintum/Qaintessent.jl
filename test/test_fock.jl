using Test
using TestSetExtensions
using LinearAlgebra
using Qaintessent
# using QuantumOptics

# @testset ExtendedTestSet "fock state photon number as integer" begin
#     modes = 4
#     trunc = 5
#     s = 3
#     fs = FockState(s, modes, trunc)
#
#     fs_ref = FockState{modes,trunc}(ComplexF64[
#                                 0 0 1.0 0 0;
#                                 0 0 1.0 0 0;
#                                 0 0 1.0 0 0;
#                                 0 0 1.0 0 0], true)
#     @test fs ≈ fs_ref
# end
#
# @testset ExtendedTestSet "fock state photon number as integer exceptions" begin
#     modes = 4
#     trunc = 5
#     s = 2
#
#     neg_modes = -1
#     @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, neg_modes, trunc)
#
#     neg_trunc = -1
#     @test_throws ErrorException("Truncated photon number must be positive") FockState(s, modes, neg_trunc)
#
#     s_neg = -5
#     @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)
#
#     s_large = 7
#     @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)
# end
#
# @testset ExtendedTestSet "fock state photon number as array" begin
#     modes = 4
#     trunc = 5
#     s = [1,2,3,4]
#     fs = FockState(s, modes, trunc)
#
#     fs_ref = FockState{modes,trunc}(ComplexF64[
#                                 1.0 0 0 0 0;
#                                 0 1.0 0 0 0;
#                                 0 0 1.0 0 0;
#                                 0 0 0 1.0 0], true)
#     @test fs ≈ fs_ref
# end
#
# @testset ExtendedTestSet "fock state exceptions" begin
#     modes = 4
#     trunc = 5
#     s = [1,2,3,4]
#
#     neg_modes = -1
#     @test_throws ErrorException("Number of modes in Fock state must be positive") FockState(s, neg_modes, trunc)
#
#     neg_trunc = -1
#     @test_throws ErrorException("Truncated photon number must be positive") FockState(s, modes, neg_trunc)
#
#     s_short = [1,2,4]
#     @test_throws ErrorException("4-mode FockState created but 3 modes provided in `s`") FockState(s_short, modes, trunc)
#
#     s_neg = [-1,2,3,4]
#     @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)
#
#     s_large = [8,2,3,4]
#     @test_throws ErrorException("Photon number must be positive and less than truncated photon number") FockState(s_neg, modes, trunc)
# end
#
# @testset ExtendedTestSet "annihilation and creation operators" begin
#     modes = 3
#     trunc = 4
#     s = [1,2,3]
#     fs = FockState(s, modes, trunc)
#
#     cg = â{trunc}()
#     fs2 = apply(cg, deepcopy(fs))
#     fs_ref = FockState{modes,trunc}(ComplexF64[
#                                 0 0 0 0;
#                                 1.0 0 0 0;
#                                 0 sqrt(2) 0 0], true)
#     @test fs2 ≈ fs_ref
#
#     cg = âDag{trunc}()
#     fs2 = apply(cg, deepcopy(fs))
#     fs_ref = FockState{modes,trunc}(ComplexF64[
#                                 0 1.0 0 0;
#                                 0 0 sqrt(2) 0;
#                                 0 0 0 sqrt(3)], true)
#     @test fs2 ≈ fs_ref
# end
#
# @testset ExtendedTestSet "photon number operator" begin
#     modes = 3
#     trunc = 4
#     s = [4,1,3]
#     fs = FockState(s, modes, trunc)
#
#     cg = n̂{trunc}()
#     fs2 = apply(cg, deepcopy(fs))
#     fs_ref = FockState{modes,trunc}(ComplexF64[
#                                 0 0 0 4;
#                                 1 0 0 0;
#                                 0 0 3 0], true)
#     @test fs2 ≈ fs_ref
# end

# @testset ExtendedTestSet "coherent states" begin
#     modes = 3
#     trunc = 8
#
#     # alpha = (0.3+0.2im)/norm(0.3+0.2im)
#     alpha = randn(ComplexF64, 1)[1]
#     # alpha = alpha/norm(alpha)
#
#     cg1 = D̂{trunc}(alpha)
#     cg3 = D̂{trunc}(-alpha)
#     cg2 = â{trunc}()
#
#     # println(norm(D1))
#     D = Qaintessent.matrix(cg1)
#     Ddag = Qaintessent.matrix(cg3)
#     # println(norm(D2))
#     a = Qaintessent.matrix(cg2)
#
#     ψ = zeros(ComplexF64, (trunc,))
#     ψ[1] = 1
#     ψ = D*ψ
#     println([norm(a)^2 for a in ψ])
#     println(sum([norm(a)^2 for a in ψ]))
#     ψ2 = a*ψ
#
#     α = ψ2./ψ
#     β = α[1]
#     println(β)
#     α = α/β
#     println([norm(a) for a in α])
# end
# #
# @testset ExtendedTestSet "phase-shift operator" begin
#     modes = 3
#     trunc = 5
#
#     theta = 0.34
#
#     cg1 = D̂{trunc}(theta)
#     cg2 = â{trunc}()
#
#     D1 = adjoint(Qaintessent.matrix(cg1))
#     # println(norm(D1))
#     D2 = Qaintessent.matrix(cg1)
#     # println(norm(D2))
#     a = Qaintessent.matrix(cg2)
#
#     ψ = zeros(ComplexF64, (trunc,))
#     ψ[1] = 1
#     ψ = D2*ψ
#     ψ2 = [abs(a^2) for a in ψ]
#     println(ψ2)
#     println(norm(D2*ψ))
# end
#
# @testset ExtendedTestSet "displacement operator" begin
#     modes = 3
#     trunc = 5
#
#     alpha = 1im
#
#     cg1 = D̂{trunc}(alpha)
#     cg2 = â{trunc}()
#
#     D1 = adjoint(Qaintessent.matrix(cg1))
#     # println(norm(D1))
#     D2 = Qaintessent.matrix(cg1)
#     # println(norm(D2))
#     a = Qaintessent.matrix(cg2)
#
#     ψ = [1, 0, 0, 0, 0]
#     ψ = D2*ψ
#     ψ2 = [abs(a^2) for a in ψ]
#     println(ψ2)
#     println(norm(D2*ψ))
# end
# 
# @testset ExtendedTestSet "single mode squeezing operator" begin
#     modes = 3
#     trunc = 5
#
#     alpha = 1im
#
#     cg1 = D̂{trunc}(alpha)
#     cg2 = â{trunc}()
#
#     D1 = adjoint(Qaintessent.matrix(cg1))
#     # println(norm(D1))
#     D2 = Qaintessent.matrix(cg1)
#     # println(norm(D2))
#     a = Qaintessent.matrix(cg2)
#
#     ψ = [1, 0, 0, 0, 0]
#     ψ = D2*ψ
#     ψ2 = [abs(a^2) for a in ψ]
#     println(ψ2)
#     println(norm(D2*ψ))
# end
#
# @testset ExtendedTestSet "two mode squeezing operator" begin
#     modes = 3
#     trunc = 5
#
#     alpha = 1im
#
#     cg1 = D̂{trunc}(alpha)
#     cg2 = â{trunc}()
#
#     D1 = adjoint(Qaintessent.matrix(cg1))
#     # println(norm(D1))
#     D2 = Qaintessent.matrix(cg1)
#     # println(norm(D2))
#     a = Qaintessent.matrix(cg2)
#
#     ψ = [1, 0, 0, 0, 0]
#     ψ = D2*ψ
#     ψ2 = [abs(a^2) for a in ψ]
#     println(ψ2)
#     println(norm(D2*ψ))
# end
