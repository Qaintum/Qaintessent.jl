"""
    ZXNode{M}

Construct a Vertex in ZX calculus. Contains in and out vertices
in::Vector
    Incoming legs

out::Vector
    Outoing legs
"""
mutable struct ZXNode{M}
    angle::Float64
    in::Vector{Ref{ZXNode}}
    out::Vector{Ref{ZXNode}}
end

Base.Ref(n::ZXNode{M} where {M}) = Base.Ref{ZXNode}(n)

function zx_decompose(cg::CircuitGate{M,N,HadamardGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{HadamardGate}(0, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,XGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{XGate}(0, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,YGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    wires = ZXNode[ZXNode{XGate}(0, input, Ref{ZXNode}[]), ZXNode{YGate}(0, Ref{ZXNode}[], Ref{ZXNode}[])]

    for i in 1:1
        push!(wires[i+1].in, Ref(wires[i]))
        push!(wires[i].out, Ref(wires[i+1]))
    end
    wires
end

function zx_decompose(cg::CircuitGate{M,N,ZGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(0, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,RxGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{XGate}(cg.gate.θ[1], input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,RyGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    wires = ZXNode[ZXNode{ZGate}(-π/2, input, Ref{ZXNode}[]), ZXNode{XGate}(-cg.gate.θ[1], Ref{ZXNode}[], Ref{ZXNode}[]), ZXNode{ZGate}(π/2, Ref{ZXNode}[], Ref{ZXNode}[])]
    for i in 1:2
        push!(wires[i+1].in, Ref(wires[i]))
        push!(wires[i].out, Ref(wires[i+1]))
    end
    wires
end

function zx_decompose(cg::CircuitGate{M,N,RzGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(cg.gate.θ[1], input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,RotationGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    θ1 = atan(cg.gate.n[1]/cg.gate.n[2])
    θ2 = atan((cg.gate.n[1]*sin(θ1)+cg.gate.n[2]*cos(θ1))/cg.gate.n[3])
    wires = ZXNode[ZXNode{ZGate}(θ1, input, Ref{ZXNode}[]), ZXNode{XGate}(θ2, Ref{ZXNode}[], Ref{ZXNode}[]), ZXNode{ZGate}(cg.gate.θ[], Ref{ZXNode}[], Ref{ZXNode}[]), ZXNode{XGate}(-θ2, Ref{ZXNode}[], Ref{ZXNode}[]), ZXNode{ZGate}(-θ1, Ref{ZXNode}[], Ref{ZXNode}[])]
    for i in 1:4
        push!(wires[i+1].in, Ref(wires[i]))
        push!(wires[i].out, Ref(wires[i+1]))
    end
end

function zx_decompose(cg::CircuitGate{M,N,TGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(π/4, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,TdagGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(-π/4, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,SGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(π/2, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,SdagGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(-π/2, input, Ref{ZXNode}[])]
end

function zx_decompose(cg::CircuitGate{M,N,ControlledGate}, input=Ref{ZXNode}[]) where {M,N}
    @assert length(input) ≤ M
    [ZXNode{ZGate}(-π/2, input, Ref{ZXNode}[])]
end

mutable struct ZXCircuit
    factor::Float64
    vertices::Vector{Ref{ZXNode}}

    @doc """
        ZXCircuit(cgc::CircuitGateChain)
    Generates ZXCircuit object from CircuitGateChain.
    """
    function ZXCircuit(cgc::CircuitGateChain{N}) where {N}
        wires = Vector{Ref{ZXNode}}(undef, N)
        nodes = Ref{ZXNodes}[]
        for moment in cgc
            for gate in moment
                input = copy(wires[gate.iwires])
                gates = zx_decompose!(gate, input)
                append!(nodes, gates)
                wires[gate.iwires] = output(gates)
            end
        end
        new(1, nodes)
    end

    @doc """
        ZXCircuit(cg::CircuitGate)
    Generates ZXCircuit object from CircuitGate.
    """
    function ZXCircuit(cg::CircuitGate)
        new(0, zx_decompose(cg))
    end
end

function Base.getindex(zx::ZXCircuit, i::Integer)
    1 <= i <= length(zx.vertices) || throw(BoundsError(zx, i))
    return zx.vertices[i]
end

function Base.iterate(zx::ZXCircuit, state=1) where {N}
    return state > length(zx.vertices) ? nothing : (zx[state], state+1)
end

# implement methods required for iteration
function Base.firstindex(zx::ZXCircuit)
    return 1
end

function Base.lastindex(zx::ZXCircuit)
    return length(zx.vertices)
end

function Base.length(zx::ZXCircuit)
    return length(zx.vertices)
end
