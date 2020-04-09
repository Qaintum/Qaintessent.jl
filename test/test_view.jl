using Test
using TestSetExtensions
using LinearAlgebra
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
    printed_output = String[]
    correct_output = ["",
                      "    1 ————————•—————x————[Ry]———x———",
                      "              |     |           |   ",
                      "    2 ———————[Rx]———•————[Pϕ]———x———",
                      "              |     |               ",
                      "    3 —[H ]——————————————[Rϕ]———x———",
                      "              |     |           |   ",
                      "    4 ————————•—————•———————————•———",
                      "                    |           |   ",
                      "    5 ——————————————x———————————x———"]

    original_stdout = stdout
    (rd, wr) = redirect_stdout()
    print(cgc)
    for line in correct_output
        push!(printed_output, readline(rd))
    end
    redirect_stdout(original_stdout)
    @test all(isequal(correct_output, printed_output))

end
