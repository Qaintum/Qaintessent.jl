using Test
using TestSetExtensions
using Qaintessent


@testset ExtendedTestSet "test view" begin

    N = 5
    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate(2, (1, 4), RxGate(√0.2), N),
        controlled_circuit_gate((1,5), (2,4), SwapGate(), N),
        single_qubit_circuit_gate(2, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N),
        two_qubit_circuit_gate(1,2, SwapGate(), N),
        controlled_circuit_gate((3,5), 4, SwapGate(), N),
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


@testset ExtendedTestSet "test view with moments" begin

    N = 5
    cgc = CircuitGateChain{N}([
        Moment{N}(
        [single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate(2, (1, 4), RxGate(√0.2), N)]
        ),
        Moment{N}(controlled_circuit_gate((1,5), (2,4), SwapGate(), N)),
        Moment{N}([single_qubit_circuit_gate(2, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N)]),
        Moment{N}([two_qubit_circuit_gate(1,2, SwapGate(), N),
        controlled_circuit_gate((3,5), 4, SwapGate(), N)]),
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

@testset ExtendedTestSet "test view moments" begin

    N = 5
    m = Moment{N}([single_qubit_circuit_gate(2, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N)])

    m_refstring =
        "\n" *
        "    5 ...|——————|...\n" *
        "         |      |   \n" *
        "    4 ...|——————|...\n" *
        "         |      |   \n" *
        "    3 ...|—[Rθ]—|...\n" *
        "         |      |   \n" *
        "    2 ...|—[Pϕ]—|...\n" *
        "         |      |   \n" *
        "    1 ...|—[Ry]—|...\n" *
        "\n"
    io = IOBuffer()
    show(io, [m])
    @test String(take!(io)) == m_refstring

    m = Moment{N}([single_qubit_circuit_gate(5, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N),
        controlled_circuit_gate((4), (2), HadamardGate(), N)])

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

    m = [Moment{N}([single_qubit_circuit_gate(5, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N),
        controlled_circuit_gate((4), (2), HadamardGate(), N)]),
        Moment{N}([single_qubit_circuit_gate(5, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N),
        controlled_circuit_gate((4), (2), HadamardGate(), N)]),
        ]

    m_refstring =
        "Array{Moment{5},1}[\n" *
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

    m = Moment{N}([single_qubit_circuit_gate(5, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        controlled_circuit_gate((2), (4), HadamardGate(), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N)])
    m_refstring =
    "CircuitGate{1,5,PhaseShiftGate}((5,), PhaseShiftGate([0.6283185307179586]), Int64[])\n" *
    "CircuitGate{1,5,RotationGate}((3,), RotationGate([0.3141592653589793, 0.0, 0.0]), Int64[])\n" *
    "CircuitGate{2,5,ControlledGate{1,2,HadamardGate}}((2, 4), ControlledGate{1,2,HadamardGate}(HadamardGate()), Union{Int64, Expr}[])\n" *
    "CircuitGate{1,5,RyGate}((1,), RyGate([4.39822971502571]), Int64[])\n"

    io = IOBuffer()
    show(io, m)
    @test String(take!(io)) == m_refstring
end
