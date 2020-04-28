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
    ["—[Rθ]—"], i
end

function view(g::PhaseShiftGate, i::Vector{Int64})
    ["—[Pϕ]—"], i
end

function view(g::SGate, i::Vector{Int64})
    ["—[S ]—"], i
end

function view(g::TGate, i::Vector{Int64})
    ["—[T ]—"], i
end

function view(g::SdagGate, i::Vector{Int64})
    ["—[S†]—"], i
end

function view(g::TdagGate, i::Vector{Int64})
    ["—[T†]—"], i
end

function view(g::SwapGate, i::Vector{Int64})
    ["——x———", "——x———"], i
end

view(g) = view(g, [1])

Base.size(g::AbstractGate{N}) where {N} = N

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

function updatecol(g::CircuitGate{M,N,G}, nw::AbstractVector{String}, ng::AbstractVector{String}) where {M,N,G}
    iwire = [x for x in g.iwire]
    if length(g.iwire) > 1
        index = min(iwire...):max(iwire...)-1
        ng[index] .= "  |   "
    end
    gate, index = view(g.gate, iwire)
    nw[index] .= gate
    return nw, ng
end

function wire_enum(N::Int)
    wires = lpad.(1:N, 5, " ")
    wires = rpad.(wires, 6, " ")
    gaps = fill("      ", N-1)
    return wires, gaps
end

function Base.show(io::IO, m::Moment{N}) where {N}
    w, g = wire_enum(N)
    i = Int[]
    w = w .* "...|"
    g = g .* "   |"
    nw = fill("——————", N)
    ng = fill("      ", N-1)
    for gate in m
        r = min(gate.iwire...):max(gate.iwire...)
        if length(intersect(i, r)) == 0
            i = union(i, r)
            nw, ng = updatecol(gate, nw, ng)
        else
            i = r
            w = w .* nw
            g = g .* ng
            nw = fill("——————", N)
            ng = fill("      ", N-1)
            nw, ng = updatecol(gate, nw, ng)
        end
    end
    w = w .* nw .* "|..."
    g = g .* ng .* "|   "
    circuit = interleave(w, g)
    println(io, "")
    for c in circuit
        println(io, c)
    end
end

function Base.show(io::IO, c::CircuitGateChain{N}) where {N}
    w, g = wire_enum(N)
    i = Int[]
    nw = fill("——————", N)
    ng = fill("      ", N-1)
    for moment in c
        for gate in moment
            r = min(gate.iwire...):max(gate.iwire...)
            if length(intersect(i, r)) == 0
                i = union(i, r)
                nw, ng = updatecol(gate, nw, ng)
            else
                i = r
                w = w .* nw
                g = g .* ng
                nw = fill("——————", N)
                ng = fill("      ", N-1)
                nw, ng = updatecol(gate, nw, ng)
            end
        end
    end
    w = w .* nw
    g = g .* ng
    circuit = interleave(w, g)
    println(io, "")
    for c in circuit
        println(io, c)
    end
end

Base.show(io::IO, c::Circuit) = Base.show(io, c.cgc)
