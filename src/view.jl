
function view(::HadamardGate, i::Vector{Int64})
    ["—[H ]—"], i
end

function view(::XGate, i::Vector{Int64})
    ["—[X ]—"], i
end

function view(::YGate, i::Vector{Int64})
    ["—[Y ]—"], i
end

function view(::ZGate, i::Vector{Int64})
    ["—[Z ]—"], i
end

function view(::RxGate, i::Vector{Int64})
    ["—[Rx]—"], i
end

function view(::RyGate, i::Vector{Int64})
    ["—[Ry]—"], i
end

function view(::RzGate, i::Vector{Int64})
    ["—[Rz]—"], i
end

function view(::RotationGate, i::Vector{Int64})
    ["—[Rθ]—"], i
end

function view(::PhaseShiftGate, i::Vector{Int64})
    ["—[Pϕ]—"], i
end

function view(::SGate, i::Vector{Int64})
    ["—[S ]—"], i
end

function view(::TGate, i::Vector{Int64})
    ["—[T ]—"], i
end

function view(::SdagGate, i::Vector{Int64})
    ["—[S†]—"], i
end

function view(::TdagGate, i::Vector{Int64})
    ["—[T†]—"], i
end

function view(::SwapGate, i::Vector{Int64})
    ["——x———", "——x———"], i
end

view(g) = view(g, [1])

function view(g::ControlledGate, i::Vector{Int64})
    gate = fill("——•———", length(i))
    gate[1:num_wires(g.U)] = view(g.U)[1]
    return gate, i
end

function view(::AbstractGate, i::Vector{Int64})
    fill("——□———", length(i)), i
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

function updatecol(g::CircuitGate{M,G}, nw::Vector{String}, ng::Vector{String}) where {M,G}
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

function momentdiagram(io::IO, m::Moment, N::Int)
    if N == 0
        return
    end
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
    for c in reverse(circuit)
        println(io, c)
    end
    println(io, "")
end

function Base.show(io::IO, m::Vector{Moment})
    if length(m) == 0
        return
    end
    N = maximum(req_wires.(m))
    for moment in m
        momentdiagram(io, moment, N)
    end
end

function Base.show(io::IO, m::Moment)
    for gate in m
        println(io, gate)
    end
end

function Base.show(io::IO, cgs::Vector{<:CircuitGate}, N::Union{Nothing,Int}=nothing)
    if isempty(cgs)
        println(io, "[]")
        return
    end
    if isnothing(N)
        N = maximum(req_wires.(cgs))
    end
    w, g = wire_enum(N)
    i = Int[]
    nw = fill("——————", N)
    ng = fill("      ", N-1)
    for cg in cgs
        r = min(cg.iwire...):max(cg.iwire...)
        if length(intersect(i, r)) == 0
            i = union(i, r)
            nw, ng = updatecol(cg, nw, ng)
        else
            i = r
            w = w .* nw
            g = g .* ng
            nw = fill("——————", N)
            ng = fill("      ", N-1)
            nw, ng = updatecol(cg, nw, ng)
        end
    end
    w = w .* nw
    g = g .* ng
    circuit = interleave(w, g)
    println(io, "")
    for c in reverse(circuit)
        println(io, c)
    end
end

function _show(io::IO, c::Vector{Moment}, N::Int)
    cg = getproperty.(c, :gates)
    cg = collect(Base.Iterators.flatten(cg))
    show(io::IO, cg, N)
end

"""
    Base.show(io::IO, c::Circuit) = Base.show(io, c.cgc)

extends base `show` method to draw visual representation of a `Circuit{N}` object.
"""
Base.show(io::IO, c::Circuit{N}) where {N} = _show(io, c.moments, N)
