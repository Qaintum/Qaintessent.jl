
"""
    DagGate

    Construct a Vertex in the Dag,contains gate information
    cg::Union{CircuitGate,Nothing}
        circuit gate contained at vertex

    prev::Union{DagGate,Nothing}
        previous vertex in this DAG

    next::Union{DagGate,Nothing}
        next vertex in this DAG

    connected::Union{AbstractVector{Ref{DagGate}},Nothing}
        all other vertices in this multi-wire DAG
"""

mutable struct DagGate
    cg::Union{CircuitGate,Nothing}
    prev::Union{DagGate,Nothing}
    next::Union{DagGate,Nothing}
    connected::Union{AbstractVector{Ref{Union{DagGate,Nothing}}},Nothing}
end

const RefDagGate = Ref{Union{DagGate,Nothing}}
Base.Ref(d::DagGate) = Base.Ref{Union{DagGate,Nothing}}(d)
Base.Ref(d::Nothing) = Base.Ref{Union{DagGate,Nothing}}(d)

#copy constructor
DagGate(dg::DagGate) = DagGate(dg.cg::Union{CircuitGate,Nothing}, dg.prev::Union{DagGate,Nothing}, dg.next::Union{DagGate,Nothing}, dg.connected::AbstractArray{RefDagGate})

DagGate(cg::CircuitGate) = DagGate(cg::Union{CircuitGate,Nothing}, nothing, nothing, nothing)

DagGate() = DagGate(nothing, nothing, nothing, nothing)

function DagGate(dg::DagGate, cg::CircuitGate, con::Union{AbstractVector{RefDagGate},Nothing})
    while true
        if isnothing(dg.next)
            break
        end
        dg = dg.next
    end
    d = DagGate(cg,dg,nothing,con)
    dg.next = d
    return d
end

function Base.length(d::DagGate)
    count = 1
    while !isnothing(d.next)
        count += 1
        d = d.next
    end
    return count
end


"""
    Dag

    Construct a Directed Acyclic Graph to represent a cgc
    iwire::Union{AbstractVector{DagGate},Nothing}
        vector of all wires in circuit
"""

mutable struct Dag
    iwire::AbstractVector{Union{DagGate,Nothing}}
    function Dag(size::Integer)
        dgs = DagGate[]
        for i in 1:size
            push!(dgs,DagGate())
        end
        new(dgs)
    end


    function Dag(cgc::CircuitGateChain{N}) where {N}
        dgs = DagGate[]
        for i in 1:N
            push!(dgs,DagGate())
        end

        for moment in cgc
            for gate in moment
                con = RefDagGate[]
                for wire in gate.iwire
                    d = DagGate(dgs[wire],gate,con)
                    push!(con,Ref(d))
                end
            end
        end

        new(dgs)
    end

end

function Base.size(d::Dag)
    count = 0
    for i in 1:length(d.iwire)
        wire = d.iwire[i]
        while wire.next != nothing
            count += 1
            wire = wire.next
        end
    end
    count
end

"""
    remove!
        d::DagGate

    basic dag operations to move gates around. removes d from given Dag representation

"""

function remove!(d::RefDagGate)
    !isnothing(d[].cg) || error("Cannot remove empty gate!")
    if isnothing(d[].next)
        d[].prev.next = nothing
    else
        prev = d[].prev
        d[].prev.next = d[].next
        d[].next.prev = prev
    end

    return d
end

"""
    insert!
        loc::RefDagGate
            position to insert gate after

        d::DagGate
            gate to be inserted

    basic dag operations to move gates around. inserts d after loc

"""

function insert!(loc::RefDagGate, d::RefDagGate)
    if isnothing(loc[].next)
        loc[].next = d[]
        d[].prev = loc[]
    else
        d[].prev = loc[]
        d[].next = loc[].next
        loc[].next = d[]
        d[].next.prev = d[]
    end
    return loc
end


"""
    firstinsert!
        loc::RefDagGate
            position to insert gate before

        d::DagGate
            gate to be inserted

    basic dag operations to move gates around. inserts d before loc

"""

function firstinsert!(loc::RefDagGate, d::RefDagGate)
    copy = DagGate(loc[].cg)
    insert!(loc::RefDagGate,Ref(copy))
    loc[].cg = d[].cg
    return loc
end

"""
    get_controls
        cg::CircuitGate

    get control wires from circuit gate

"""

function get_controls(cg::CircuitGate{N,M,G}) where {N,M,G <:ControlledGate{O,P}} where {O,P}
    num_gate_wires = O
    num_total_wires = P
    iwire = cg.iwire
    cntrl = iwire[1:P-O]
    gate = iwire[end-O+1:end]
    (cntrl,gate)
end

function get_controls(cg::CircuitGate)
    ((),cg.iwire)
end

"""
    get_total_wires
        cg::CircuitGate

    get total wires from circuit gate
"""

function get_total_wires(cg::CircuitGate{N,M,G}) where {N,M,G}
    return M
end

"""
    check_comm
        d1::RefDagGate
        d2::RefDagGate

    function to check if two Vertices in DAG can commute to each other. this
    methods checks to see if d1 can commute backwards to d2.
    note that d1 should be in front of d2
"""

function check_comm(d1::RefDagGate, d2::RefDagGate)
    temp = d1[].prev
    d1cg = d1[].cg
    while !(d2[].cg === temp.cg)
        if !iscommuting(temp.cg,d1cg)
            return false
        end
        temp = temp.prev
    end
    return true
end


"""
    match
        dag1::DagGate
        dag2::DagGate

    check if given DagGate sequence matches a pattern

"""

function match(g1::ControlledGate{N}, g2::ControlledGate{N}) where {N}
    if typeof(g1.U) != typeof(g2.U)
        return false
    end
    return true
end

function match(g1::AbstractGate{N}, g2::AbstractGate{N}) where {N}
    if typeof(g1) != typeof(g2)
        return false
    end
    return true
end

match(g1::AbstractGate{M}, g2::AbstractGate{N}) where {M,N} = false


function match!(dag_ref::RefDagGate,pattern_ref::DagGate)
    if !match(dag_ref[].cg.gate,pattern_ref.cg.gate)
        return false
    end

    total_wires = get_total_wires(dag_ref[].cg)
    refs = Vector{RefDagGate}(undef,total_wires)
    matched = RefDagGate[]
    for wire in dag_ref[].cg.iwire
        refs[wire] = dag_ref
    end

    if isnothing(dag_ref[].next)
        return false
    end

    dag = Ref(dag_ref[].next)
    pattern = pattern_ref.next


    while true
        if !match(dag[].cg.gate,pattern.cg.gate)
            while true
                if !isnothing(dag[].next)
                    dag = Ref(dag[].next)
                else
                    return false
                end

                if match(dag[].cg.gate,pattern.cg.gate)
                    for (i,wire) in enumerate(dag[].cg.iwire)
                        if isassigned(refs,wire)
                            check_comm(dag[].connected[i],refs[wire]) || return false
                        end
                    end
                    for (i,wire) in enumerate(dag[].cg.iwire)
                        refs[i] = dag
                    end
                    push!(matched,dag)
                    break
                end
            end
        else
            push!(matched,dag)
        end

        if isnothing(pattern.next)
            break
        elseif isnothing(dag[].next)
            return false
        end

        dag = Ref(dag[].next)
        pattern = pattern.next
    end

    @assert length(matched) == length(pattern_ref)-1
    for dag in matched
        dwire = collect(dag[].cg.iwire)
        refwire = collect(dag_ref[].cg.iwire)
        intersected_wires = intersect(dwire,refwire)
        for wire in intersected_wires
            d1 = dag[].connected[indexin([wire],dwire)[1]]
            d2 = dag_ref[].connected[indexin([wire],refwire)[1]]
            if !(d1[].prev.cg === d2[].cg)
                insert!(d2,remove!(d1))
            end
        end
        dag_ref = Ref(dag_ref[].next)
    end
    return true
end


"""
    opt_hadamard
        dag::Dag

    optimizes Dag object by attempting to remove Hadamard gates. based on optimization algorithm from arXiv:1710.07345v2

"""

"""
N = 1
"""

hadamard_inverse_cgc = CircuitGateChain{1}([single_qubit_circuit_gate(1, HadamardGate(), 1),
                                single_qubit_circuit_gate(1, HadamardGate(), 1)])
hadamard_inverse = Dag(hadamard_inverse_cgc).iwire[1].next

hxh_cgc = CircuitGateChain{1}([single_qubit_circuit_gate(1, HadamardGate(), 1),
                                single_qubit_circuit_gate(1, X, 1),
                                single_qubit_circuit_gate(1, HadamardGate(), 1)])
hxh = Dag(hxh_cgc).iwire[1].next
hzh_cgc = CircuitGateChain{1}([single_qubit_circuit_gate(1, HadamardGate(), 1),
                                single_qubit_circuit_gate(1, Z, 1),
                                single_qubit_circuit_gate(1, HadamardGate(), 1)])
hzh = Dag(hzh_cgc).iwire[1].next
"""
N = 2
"""
hcxh_cgc = CircuitGateChain{2}([single_qubit_circuit_gate(1, HadamardGate(), 2),
                                controlled_circuit_gate((2), 1, X, 2),
                                single_qubit_circuit_gate(1, HadamardGate(), 2)])
hcxh = Dag(hcxh_cgc).iwire[1].next

hczh_cgc = CircuitGateChain{2}([single_qubit_circuit_gate(1, HadamardGate(), 2),
                                controlled_circuit_gate((2), 1, Z, 2),
                                single_qubit_circuit_gate(1, HadamardGate(), 2)])
hczh = Dag(hczh_cgc).iwire[1].next

function hadamard_inverse_opt(d::RefDagGate)
    d = remove!(remove!(d))
end

function hxh_opt(d::RefDagGate)
    _,wires = get_controls(d[].cg)
    cg = single_qubit_circuit_gate(wires[1],Z,length(wires))
    d[].cg = cg
    remove!(remove!(Ref(d[].next)))
    return d
end

function hzh_opt(d::RefDagGate)
    _,wires = get_controls(d[].cg)
    cg = single_qubit_circuit_gate(wires[1],X,length(wires))
    d[].cg = cg
    remove!(remove!(Ref(d[].next)))
    return d
end

function hcxh_opt(d::RefDagGate)
    cntrl,wires = get_controls(d[].next.cg)

    con = RefDagGate[]
    cg = controlled_circuit_gate(cntrl[1],wires[1],Z,length(d[].next.cg.iwire))

    d[].next.cg = cg
    c = d[].next.connected[1]
    c[].cg = cg

    remove!(Ref(d[].next.next))
    remove!(d)
end

function hczh_opt(d::RefDagGate)
    cntrl,wires = get_controls(d[].next.cg)

    con = RefDagGate[]
    cg = controlled_circuit_gate(cntrl[1],wires[1],X,length(d[].next.cg.iwire))

    d[].next.cg = cg
    c = d[].next.connected[1]
    c[].cg = cg

    remove!(Ref(d[].next.next))
    remove!(d)
end

hadamard_patterns = IdDict(hadamard_inverse => hadamard_inverse_opt,
                            hxh => hxh_opt,
                            hzh => hzh_opt,
                            hcxh => hcxh_opt,
                            hczh => hczh_opt,
                        )

function opt_hadamard(dag::Dag)
    for i in 1:length(dag.iwire)
        d = Ref(Ref(dag.iwire,i)[].next)
        while !isnothing(d[]) && !isnothing(d[].cg)
            if typeof(d[].cg.gate) == HadamardGate
                for pattern in keys(hadamard_patterns)
                    if match!(d,pattern)
                        hadamard_patterns[pattern](d)
                        break
                    end
                end
            end
            d = Ref(d[].next)
        end
    end
    return dag
end

"""
    opt_adjoint
        dag::Dag

    basic optimization to remove adjoint gates. based on optimization algorithm from arXiv:1710.07345v2.

"""
function opt_adjoint(dag::Dag)
    for i in 1:length(dag.iwire)
        d = Ref(Ref(dag.iwire,i)[].next)
        while !isnothing(d[]) && !isnothing(d[].cg)
            ref = Ref(d[].next)
            adj = typeof(adjoint(d[].cg))
            while !isnothing(ref[]) && !isnothing(ref[].cg)
                iscomm = true
                if typeof(ref[].cg) == adj && get_controls(ref[].cg) == get_controls(d[].cg)
                    for (rtemp,dtemp) in zip(ref[].connected,d[].connected)
                        if !check_comm(rtemp,dtemp)
                            iscomm = false
                            break
                        end
                    end
                    iscomm || break
                    for (rtemp,dtemp) in zip(ref[].connected,d[].connected)
                        remove!(rtemp)
                        remove!(dtemp)
                    end
                    break
                end
                ref = Ref(ref[].next)
            end
            d = Ref(d[].next)
        end
    end
    return dag
end


function Base.show(io::IO,daggate::DagGate)

    print("[")

    if !isnothing(daggate.cg)
        print(string(daggate.cg) * ",")
    end

    while true
        if isnothing(daggate.next)
            break
        end
        daggate = daggate.next
        if isnothing(daggate.cg)
            break
        end
        print(string(daggate.cg) * ",")
    end
    print("]")
end

function Base.show(io::IO,wires::AbstractVector{DagGate})
    for daggate in wires
        println(daggate)
    end
end

function Base.show(io::IO, dag::Dag)
    for daggate in dag.iwire
        println(daggate)
    end
end


"""
    append_gate!
        cgs::AbstractVector{CircuitGate}
        dg::DagGate
"""
    function append_gate!(cgs::AbstractVector{AbstractCircuitGate{N}}, dag::Dag,dg::DagGate) where {N}
        for j in dg.cg.iwire
            while !(dag.iwire[j].cg === dg.cg)
                dag = append_gate!(cgs,dag,dag.iwire[j])
            end
            dag.iwire[j] = dag.iwire[j].next
        end
        push!(cgs,dg.cg)
        return dag
    end

"""
    CircuitGateChain
        dag::Dag

    method to convert DAG back to a CGC
"""

function CircuitGateChain(dag_ref::Dag)
    dag = deepcopy(dag_ref)
    N = length(dag.iwire)
    cgs = AbstractCircuitGate{N}[]

    for i in 1:N
        dag.iwire[i] = dag.iwire[i].next
    end

    for i in 1:N
        while !isnothing(dag.iwire[i])
            dag = append_gate!(cgs,dag,dag.iwire[i])
        end
    end
    CircuitGateChain{N}(cgs)
end
