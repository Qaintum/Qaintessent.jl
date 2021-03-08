using Test
using TestSetExtensions
using Qaintessent


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "test view" begin

    N = 5
    cgc = [
        circuit_gate(3, HadamardGate()),
        circuit_gate(2, RxGate(√0.2), (1, 4)),
        circuit_gate((1,5), SwapGate(), (2,4)),
        circuit_gate(2, TGate()),
        circuit_gate(2, PhaseShiftGate(0.2π)),
        circuit_gate(2, SdagGate()),
        circuit_gate(3, RotationGate(0.1π, [1, 0, 0])),
        circuit_gate(1, RyGate(1.4π)),
        circuit_gate(3, RzGate(0.4π)),
        circuit_gate(1,2, SwapGate()),
        circuit_gate((3,5), SwapGate(), 4),
        circuit_gate((3,), Y, 2),
        circuit_gate((5,), Z, 1),
        circuit_gate((4,), X)
    ]

    cgc_refstring =
        "\n" *
        "    5 ——————————————x—————————————————————————————x——————————[Z ]———————\n" *
        "                    |                             |           |         \n" *
        "    4 ————————•—————•—————————————————————————————•————————————————[X ]—\n" *
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


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "test view with circuits" begin

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


##==----------------------------------------------------------------------------------------------------------------------


@testset ExtendedTestSet "test view moments" begin

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
        "Array{Moment,1}[\n" *
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
    m_refstring =
    "CircuitGate{1,PhaseShiftGate}((5,), PhaseShiftGate([0.6283185307179586]))\n" *
    "CircuitGate{1,RotationGate}((3,), RotationGate([0.3141592653589793, 0.0, 0.0]))\n" *
    "CircuitGate{2,ControlledGate{HadamardGate}}((2, 4), ControlledGate{HadamardGate}(HadamardGate(), 1))\n" *
    "CircuitGate{1,RyGate}((1,), RyGate([4.39822971502571]))\n"

    io = IOBuffer()
    show(io, m)
    @test String(take!(io)) == m_refstring
end
