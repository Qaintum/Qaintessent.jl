@info("Ensuring example environment instantiated...")
import Pkg
Pkg.activate("../")
Pkg.instantiate()

using Parameters
@info("Loading Qaintessent, Zygote...")
using Qaintessent, Flux, Zygote, LinearAlgebra
using Zygote: @adjoint

function ngradient(f, xs::AbstractArray...)
    grads = zero.(xs)
    for (x, Δ) in zip(xs, grads), i in 1:length(x)
        δ = sqrt(eps())
        tmp = x[i]
        x[i] = tmp - δ/2
        y1 = f(xs...)
        x[i] = tmp + δ/2
        y2 = f(xs...)
        x[i] = tmp
        Δ[i] = (y2-y1)/δ
    end
    return grads
end

function cgadjoint(cg, ψ, Δ)
    adjoint_dict = Qaintessent.backward(cg, ψ, Δ)
    gate = cg.gate
    if gate isa ControlledGate
        gate = cg.gate.U
    end
    fields = fieldnames(typeof(gate))
    adjoint_dict = [adjoint_dict[getfield(gate, name)] for name in fields]
    return adjoint_dict
end

function ψadjoint(cg, Δ)
    Qaintessent.apply(adjoint(cg), Δ)
end

@adjoint Qaintessent.apply(cg::CircuitGate, ψ::AbstractVector) = Qaintessent.apply(cg::CircuitGate, ψ::AbstractVector), Δ -> (cgadjoint(cg, ψ, Δ), (ψadjoint(cg, Δ)))

@adjoint function Qaintessent.ControlledGate{M,N}(U) where {M,N}
     Qaintessent.ControlledGate{M,N}(U), Δ -> (Tuple(real(sum(Δ))))
end
@adjoint function Qaintessent.CircuitGate{M,N}(iwire::NTuple, gate::Qaintessent.AbstractGate) where {M,N}
    return Qaintessent.CircuitGate{M,N}(iwire::NTuple, gate::Qaintessent.AbstractGate), Δ -> (nothing, Δ)
end

@adjoint single_qubit_circuit_gate(iwire, gate, N) = single_qubit_circuit_gate(iwire, gate, N), Δ -> (nothing, Δ, nothing)
@adjoint two_qubit_circuit_gate(iwire1, iwire2, gate, N) = two_qubit_circuit_gate(iwire1, iwire2, gate, N), Δ -> (nothing, nothing, Δ, nothing)
@adjoint controlled_circuit_gate(icntrl, itarget, U, N) = controlled_circuit_gate(icntrl, itarget, U, N), Δ -> (nothing, nothing, Δ, nothing)

# construct parametrized circuit

N = 4
meas = MeasurementOps{N}([Matrix{Float64}(I, 2^N, 2^N), Hermitian(randn(ComplexF64, 2^N, 2^N))])
# input quantum state
ψ = randn(ComplexF64, 2^N)
ψ = ψ/norm(ψ)
# fictitious gradients of cost function with respect to circuit output
Δ = [0.3, -1.2]

rz = RzGate(1.5π)
ps = PhaseShiftGate(0.3)
ry = RyGate(√2)
n = randn(Float64, 3)
n = n/norm(n)
rg = RotationGate(0.2π, n)

# Zygote assumes input gradients are [1, 1]
gs = gradient(Params([rz, ps, rg, ry])) do
    cgc = CircuitGateChain{N}([
        single_qubit_circuit_gate(3, HadamardGate(), N),
        controlled_circuit_gate((1, 4), 2, rz, N),
        two_qubit_circuit_gate(2, 3, SwapGate(), N),
        single_qubit_circuit_gate(3, ps, N),
        single_qubit_circuit_gate(3, rg, N),
        single_qubit_circuit_gate(1, ry, N),
    ])
    c = Circuit(cgc, meas)
    0.5 * dot(Δ, apply(c, ψ))
end

cgc = CircuitGateChain{N}([
    single_qubit_circuit_gate(3, HadamardGate(), N),
    controlled_circuit_gate((1, 4), 2, rz, N),
    two_qubit_circuit_gate(2, 3, SwapGate(), N),
    single_qubit_circuit_gate(3, ps, N),
    single_qubit_circuit_gate(3, rg, N),
    single_qubit_circuit_gate(1, ry, N),
])

grads = Qaintessent.gradients(Circuit(cgc, meas), ψ, Δ)

f(args...) = dot(Δ, apply(Circuit(cgc, meas), ψ))
ngrad = ngradient(f, rz.θ, ps.ϕ, rg.nθ, ry.θ)

println("Output from gradient calc: " * string(grads[rz.θ]) * " " * string(grads[ps.ϕ]) * " " * string(grads[rg.nθ]) * " " * string(grads[ry.θ]))
println("Output from ngradient calc: " * string(ngrad))
println("Output from Zygote calc: " * string(gs[rz]) * " " * string(gs[ps]) * " " * string(gs[rg]) * " " * string(gs[ry]))
