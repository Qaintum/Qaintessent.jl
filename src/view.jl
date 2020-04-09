function view(g::HadamardGate, i::Vector{Int64})
    ["—[H ]—"], i
end

function view(g::XGate, i::Vector{Int64})
    ["—[X ]—"], i
end

function view(g::YGate, i::Vector{Int64})
    ["—[Y ]—"], i
end

function view(g::ZGate, i::Vector{Int64})
    ["—[Z ]—"], i
end

function view(g::RxGate, i::Vector{Int64})
    ["—[Rx]—"], i
end

function view(g::RyGate, i::Vector{Int64})
    ["—[Ry]—"], i
end

function view(g::RzGate, i::Vector{Int64})
    ["—[Rz]—"], i
end

function view(g::RotationGate, i::Vector{Int64})
    ["—[Rϕ]—"], i
end

function view(g::PhaseShiftGate, i::Vector{Int64})
    ["—[Pϕ]—"], i
end

function view(g::TGate, i::Vector{Int64})
    ["—[T ]—"], i
end

function view(g::SGate, i::Vector{Int64})
    ["—[S ]—"], i
end

function view(g::SwapGate, i::Vector{Int64})
    ["——x———", "——x———"], i
end

view(g) = view(g, [1])

function Base.size(g::AbstractGate{N}) where {N}
    N
end

function view(g::ControlledGate, i::Vector{Int64})
    gate = fill("——•———", length(i))
    gate[end-Base.size(g.U)+1:end] = view(g.U)[1]
    return gate, i
end

function interleave(a::Vector{T}, b::Vector{T}) where {T}
    length(a)-1 == length(b) || error("Vector b must be one element shorter than Vector a")
    c = Vector{T}(undef, length(a) + length(b))
    c[1] = a[1]
    i = 1
    for (x, y) in zip(b, a[2:end])
        c[i += 1] = x
        c[i += 1] = y
    end
    return c
end

function getcol(g::CircuitGate{M, N}, w::Int) where {M,N}
    circuit =
    wires = fill("——————", w)
    gaps = fill("      ", w-1)

    iwire = [x for x in g.iwire]

    if length(g.iwire) > 1
        index = min(iwire...):max(iwire...)-1
        gaps[index] .= "  |   "
    end

    gate, index = view(g.gate, iwire)
    wires[index] .= gate

    return interleave(wires, gaps)
end

function wire_enum(N::Int)
    wires = lpad.(1:N, 5, " ")
    wires = rpad.(wires, 6, " ")
    gaps = fill("      ", N-1)
    return wires, gaps
end

function Base.print(c::CircuitGateChain{N}) where {N}
    wires, gaps = wire_enum(N)
    circuit = interleave(wires, gaps)
    for gate in c
        circuit = circuit .* getcol(gate, N)
    end
    println("")
    println.(circuit)
end

Base.println(c::CircuitGateChain) = Base.print(c)

Base.show(c::Circuit) = show(c.cgc)
