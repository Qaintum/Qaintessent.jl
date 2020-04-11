using Test
using TestSetExtensions
using Qaintessent


@testset ExtendedTestSet "test view" begin

    N = 5
    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate((1, 4), 2, RxGate(√0.2), N),
        controlled_circuit_gate((2,4), (1,5), SwapGate(), N),
        single_qubit_circuit_gate(2, PhaseShiftGate(0.2π), N),
        single_qubit_circuit_gate(3, RotationGate(0.1π, [1, 0, 0]), N),
        single_qubit_circuit_gate(1, RyGate(1.4π), N),
        two_qubit_circuit_gate(1,2, SwapGate(), N),
        controlled_circuit_gate(4, (3,5), SwapGate(), N),
    ])

    cgc_refstring =
        "\n" *
        "    1 ————————•—————x————[Ry]———x———\n" *
        "              |     |           |   \n" *
        "    2 ———————[Rx]———•————[Pϕ]———x———\n" *
        "              |     |               \n" *
        "    3 —[H ]——————————————[Rθ]———x———\n" *
        "              |     |           |   \n" *
        "    4 ————————•—————•———————————•———\n" *
        "                    |           |   \n" *
        "    5 ——————————————x———————————x———\n"

    io = IOBuffer()
    show(io, cgc)
    @test String(take!(io)) == cgc_refstring

end
