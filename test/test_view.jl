using Test
using TestSetExtensions
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------

@testset ExtendedTestSet "views" begin
    @testset "view basic" begin

        N = 5
        cgc = [
            circuit_gate(3, HadamardGate()),
            circuit_gate(2, RxGate(√0.2), (1, 4)),
            circuit_gate((1,5), SwapGate(), (2,4)),
            circuit_gate(2, TGate()),
            circuit_gate(2, PhaseShiftGate(0.2π)),
            circuit_gate(2, SdagGate()),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate(4, TdagGate()),
            circuit_gate(1, RyGate(1.4π)),
            circuit_gate(3, RzGate(0.4π)),
            circuit_gate(1,2, SwapGate()),
            circuit_gate(4, SGate()),
            circuit_gate((3,5), SwapGate(), 4),
            circuit_gate((3,), Y, 2),
            circuit_gate((5,), Z, 1),
            circuit_gate((4,), X)
        ]

        cgc_refstring =
            "\n" *
            "    5 ——————————————x—————————————————————————————x——————————[Z ]———————\n" *
            "                    |                             |           |         \n" *
            "    4 ————————•—————•————————————————[T†]——[S ]———•————————————————[X ]—\n" *
            "              |     |                             |           |         \n" *
            "    3 —[H ]——————————————————————————[Rθ]——[Rz]———x————[Y ]—————————————\n" *
            "              |     |                                   |     |         \n" *
            "    2 ———————[Rx]———•————[T ]——[Pϕ]——[S†]———x———————————•———————————————\n" *
            "              |     |                       |                 |         \n" *
            "    1 ————————•—————x————————————————[Ry]———x—————————————————•—————————\n"

        io = IOBuffer()
        show(io, cgc)
        @test String(take!(io)) == cgc_refstring
        

    end

    @testset ExtendedTestSet "empty basic view" begin
        cgc = CircuitGate[]

        cgc_refstring = "[]\n"

        io = IOBuffer()
        show(io, cgc)
        @test String(take!(io)) == cgc_refstring
    end

    ##==----------------------------------------------------------------------------------------------------------------------


    @testset "view circuit" begin

        N = 5
        cgc = Circuit{N}([
            Moment(
            [circuit_gate(3, HadamardGate()),
            circuit_gate(2, RxGate(√0.2), (1, 4))]
            ),
            Moment(circuit_gate((1,5), SwapGate(), (2,4))),
            Moment([circuit_gate(2, PhaseShiftGate(0.2π)),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate(1, RyGate(1.4π))]),
            Moment([circuit_gate(1,2, SwapGate()),
            circuit_gate((3,5), SwapGate(), 4)]),
        ])

        cgc_refstring =
            "\n" *
            "    5 ——————————————x———————————x———\n" *
            "                    |           |   \n" *
            "    4 ————————•—————•———————————•———\n" *
            "              |     |           |   \n" *
            "    3 —[H ]——————————————[Rθ]———x———\n" *
            "              |     |               \n" *
            "    2 ———————[Rx]———•————[Pϕ]———x———\n" *
            "              |     |           |   \n" *
            "    1 ————————•—————x————[Ry]———x———\n"

        io = IOBuffer()
        show(io, cgc)
        @test String(take!(io)) == cgc_refstring

    end

    @testset ExtendedTestSet "view empty circuit" begin
        N = 5
        cgc = Circuit{N}()
        cgc_refstring = "[]\n"

        io = IOBuffer()
        show(io, cgc)
        @test String(take!(io)) == cgc_refstring
    end


    ##==----------------------------------------------------------------------------------------------------------------------


    @testset "view moment" begin

        N = 5
        m = Moment([circuit_gate(2, PhaseShiftGate(0.2π)),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate(1, RyGate(1.4π))])

        m_refstring =
            "\n" *
            "    3 ...|—[Rθ]—|...\n" *
            "         |      |   \n" *
            "    2 ...|—[Pϕ]—|...\n" *
            "         |      |   \n" *
            "    1 ...|—[Ry]—|...\n" *
            "\n"
        io = IOBuffer()
        show(io, [m])
        @test String(take!(io)) == m_refstring

        m = Moment([circuit_gate(5, PhaseShiftGate(0.2π)),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate(1, RyGate(1.4π)),
            circuit_gate((4), HadamardGate(), (2))])

        m_refstring =
            "\n" *
            "    5 ...|—[Pϕ]———————|...\n" *
            "         |            |   \n" *
            "    4 ...|———————[H ]—|...\n" *
            "         |        |   |   \n" *
            "    3 ...|—[Rθ]———————|...\n" *
            "         |        |   |   \n" *
            "    2 ...|————————•———|...\n" *
            "         |            |   \n" *
            "    1 ...|—[Ry]———————|...\n" *
            "\n"

        io = IOBuffer()
        show(io, [m])
        @test String(take!(io)) == m_refstring

        m = [Moment([circuit_gate(5, PhaseShiftGate(0.2π)),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate(1, RyGate(1.4π)),
            circuit_gate((4), HadamardGate(), (2))]),
            Moment([circuit_gate(5, PhaseShiftGate(0.2π)),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate(1, RyGate(1.4π)),
            circuit_gate((4), HadamardGate(), (2))]),
            ]

        m_refstring =
            "Vector{Moment}[\n" *
            "    5 ...|—[Pϕ]———————|...\n" *
            "         |            |   \n" *
            "    4 ...|———————[H ]—|...\n" *
            "         |        |   |   \n" *
            "    3 ...|—[Rθ]———————|...\n" *
            "         |        |   |   \n" *
            "    2 ...|————————•———|...\n" *
            "         |            |   \n" *
            "    1 ...|—[Ry]———————|...\n" *
            "\n" *
            "\n" *
            "    5 ...|—[Pϕ]———————|...\n" *
            "         |            |   \n" *
            "    4 ...|———————[H ]—|...\n" *
            "         |        |   |   \n" *
            "    3 ...|—[Rθ]———————|...\n" *
            "         |        |   |   \n" *
            "    2 ...|————————•———|...\n" *
            "         |            |   \n" *
            "    1 ...|—[Ry]———————|...\n" *
            "\n]"

        io = IOBuffer()
        show(io, [m])
        @test String(take!(io)) == m_refstring

        m = Moment([circuit_gate(5, PhaseShiftGate(0.2π)),
            circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
            circuit_gate((2), HadamardGate(), 4),
            circuit_gate(1, RyGate(1.4π))])
        m_refstring64 =
        "CircuitGate{1, PhaseShiftGate}((5,), PhaseShiftGate([0.6283185307179586]))\n" *
        "CircuitGate{1, RotationGate}((3,), RotationGate([0.3141592653589793, 0.0, 0.0]))\n" *
        "CircuitGate{2, ControlledGate{HadamardGate}}((2, 4), ControlledGate{HadamardGate}(HadamardGate(), 1))\n" *
        "CircuitGate{1, RyGate}((1,), RyGate([4.39822971502571]))\n"

        m_refstring32 =
        "CircuitGate{1, PhaseShiftGate}((5,), PhaseShiftGate(Float32[0.62831855]))\n" *
        "CircuitGate{1, RotationGate}((3,), RotationGate(Float32[0.31415927, 0.0, 0.0]))\n" *
        "CircuitGate{2, ControlledGate{HadamardGate}}((2, 4), ControlledGate{HadamardGate}(HadamardGate(), 1))\n" *
        "CircuitGate{1, RyGate}((1,), RyGate(Float32[4.3982296]))\n"

        if FloatQ == Float32
            m_refstring = m_refstring32
        else
            m_refstring = m_refstring64
        end

        io = IOBuffer()
        show(io, m)
        @test String(take!(io)) == m_refstring
    end

    @testset "view empty moment" begin
        N = 5
        cgc = CircuitGate[]
        m = Moment(cgc)
        m_refstring = ""

        io = IOBuffer()
        show(io, m)
        @test String(take!(io)) == m_refstring
    end


    @testset "view empty moments" begin
        N = 5
        cgc = CircuitGate[]
        m = [Moment(cgc)]
        m_refstring = ""

        io = IOBuffer()
        show(io, m)
        @test String(take!(io)) == m_refstring
    end
end