"""
    qft_circuit(N)

Construct the quantum Fourier transform circuit for `N` qubits.
"""
function qft_circuit(N)
    main = vcat([
            [j == i ?
             single_qubit_circuit_gate(i, HadamardGate(), N) :
             controlled_circuit_gate(j, i, PhaseShiftGate(2π/2^(j-i+1)), N) for j in i:N]
        for i in 1:N]...)
    swap = [two_qubit_circuit_gate(i, N-i+1, SwapGate(), N) for i in 1:(N÷2)]
    CircuitGateChain{N}(N > 1 ? [main; swap] : main)
end

"""
    qfa_circuit(N)

Construct the quantum fourier addition circuit for `N` qubits
"""
function qfa_circuit(N)
    main = vcat([[controlled_circuit_gate(j, i+N, PhaseShiftGate(1/(2^i)), 2*N) for j in i:N] for i in 1:N]...)
    CircuitGateChain{2*N}(main)
end

"""
    toffoli_circuit(cntrl, trg, N)

Construct the decomposed toffoli circuit for `N` qubits based on Nielson and Chung's decomposition
"""

function toffoli_circuit(cntrl::Tuple{Int, Int} , trg::Int, N::Int)
    main = [
            single_qubit_circuit_gate(trg, HadamardGate(),N),
            controlled_circuit_gate(cntrl[2], trg, X, N),
            single_qubit_circuit_gate(trg, TdagGate(),N),
            controlled_circuit_gate(cntrl[1], trg, X, N),
            single_qubit_circuit_gate(trg, TGate(),N),
            controlled_circuit_gate(cntrl[2], trg, X, N),
            single_qubit_circuit_gate(trg, TdagGate(),N),
            controlled_circuit_gate(cntrl[1], trg, X, N),
            single_qubit_circuit_gate(trg, TGate(),N),
            single_qubit_circuit_gate(cntrl[2], TGate(),N),
            controlled_circuit_gate(cntrl[1], cntrl[2], X, N),
            single_qubit_circuit_gate(trg, HadamardGate(),N),
            single_qubit_circuit_gate(cntrl[1], TGate(),N),
            single_qubit_circuit_gate(cntrl[2], TdagGate(),N),
            controlled_circuit_gate(cntrl[1], cntrl[2], X, N),
            ]
    CircuitGateChain{N}(main)
end

"""
    mod5_circuit(N)

Construct the mod5 circuit for `5` qubits. Based on circuit on
    D. Maslov, “Reversible logic synthesis http://webhome.cs.uvic.ca/dmaslov/“
"""
function mod5_circuit()
    N = 5
    main1 = CircuitGateChain{N}([
        controlled_circuit_gate(1, 3, X, N),
        controlled_circuit_gate(2, 4, X, N),
        ])
    main2 = toffoli_circuit((3,4), 5, N)
    main3 = CircuitGateChain{N}([
        controlled_circuit_gate(3, 4, X, N),
        controlled_circuit_gate(4, 5, X, N),
        ])

    return main1*main2*main3
end

"""
    vbe_adder_circuit(N)

Construct an in-place adder for 2 integers represented by `N` qubits.
Based on ripple-carry adder circuit on Vedral et. al in quant-ph/9511018“
    Returns a CircuitGateChain{3N+1} as there are N+1 ancillary wires
    If the two added integers are represented as:
        A = a0*1 + a1*2 + a2*4 + .. + aN*2^N
        B = b0*1 + b1*2 + b2*4 + .. + bN*2^N
    The input index should be a0a1a2..aNb0b1b2..bN+1
    i.e.
        input = fill(0, 2^(3N+1))
        input[index+1] = 1

    The output index will be in the form a0a1a2...aNc0c1c2...cN where
        C = (A+B) % (2^N + 1) = c0*1 + c1*2 + c2*4 + ... + cN*2^N
"""
function vbe_adder_circuit(N)
    M = 3*N + 1
    carry(a,b,c,d) = CircuitGateChain{M}(
    [
        controlled_circuit_gate((b, c), d, X, M),
        controlled_circuit_gate(b, c, X, M),
        controlled_circuit_gate((a, c), d, X, M),
    ])
    rcarry(a,b,c,d) = CircuitGateChain{M}(
    [
        controlled_circuit_gate((a, c), d, X, M),
        controlled_circuit_gate(b, c, X, M),
        controlled_circuit_gate((b, c), d, X, M),
    ])
    sum(a,b,c) = CircuitGateChain{M}(
    [
        controlled_circuit_gate(b, c, X, M),
        controlled_circuit_gate(a, c, X, M),
    ])
    rsum(a,b,c) = CircuitGateChain{M}(
    [
        controlled_circuit_gate(a, c, X, M),
        controlled_circuit_gate(b, c, X, M),
    ])

    cgc = carry(N+1, M, M-N, N)
    for i in 2:N
        cgc = cgc*carry(N+2-i, M+1-i, M+1-N-i, N+1-i)
    end

    cgc = cgc*CircuitGateChain{M}([controlled_circuit_gate(M-N+1, M-2N+1, X, M)]) *
            sum(2, M-N+1, M-2N+1)
    for i in N-1:-1:1
        cgc = cgc*rcarry(N+2-i, M+1-i, M+1-N-i, N+1-i) * sum(N+2-i, M+1-i, M+1-N-i)
    end
    return cgc
end


"""
    p_round(P, M, N)
        execute a P round for qcla algorithms
"""
function p_round(P::Function, M::Int, N::Int)
    # Add gates for P rounds
    pchain = CircuitGate{<:Any, M, <:Any}[]
    tmax = floor(Int, log(2, N)) - 1
    for t in 1:tmax
        mmax = floor(Int, N/2^t) - 1
        for m in 1:mmax
            push!(pchain, controlled_circuit_gate((P(t-1, 2m), P(t-1, 2m+1)), P(t, m), X, M))
        end
    end
    return pchain
end

function p_round(P::Function, M::Int, N::Int, Ñ::Int)
    # Add gates for P rounds
    pchain = CircuitGate{<:Any, M, <:Any}[]
    tmax = floor(Int, log(2, Ñ)) - 1
    for t in 1:tmax
        mmax = floor(Int, Ñ/2^t) - 1
        for m in 1:mmax
            push!(pchain, controlled_circuit_gate((P(t-1, 2m), P(t-1, 2m+1)), P(t, m), X, M))
        end
    end
    return pchain
end

"""
    g_round(P, M, N)
        execute a P round for qcla algorithms
"""
function g_round(P::Function, G::Function, M::Int, N::Int)
    # Add gates for G rounds
    gchain = CircuitGate{<:Any, M, <:Any}[]
    tmax = floor(Int, log(2, N))
    for t in 1:tmax
        mmax = floor(Int, N/2^t) - 1
        for m in 0:mmax
            push!(gchain, controlled_circuit_gate(( G(2^t*m+2^(t-1)) , P(t-1, 2m+1)), G(2^t*m + 2^t), X, M))
        end
    end
    return gchain
end

function g_round(P::Function, G::Function, M::Int, N::Int, Ñ::Int)
    # Add gates for G rounds
    gchain = CircuitGate{<:Any, M, <:Any}[]
    tmax = floor(Int, log(2, Ñ))
    for t in 1:tmax
        mmax = floor(Int, Ñ/2^t) - 1
        for m in 0:mmax
            push!(gchain, controlled_circuit_gate(( G(2^t*m+2^(t-1)) , P(t-1, 2m+1)), G(2^t*m + 2^t), X, M))
        end
    end
    return gchain
end

"""
    c_round(P, G, M, N)
"""
function c_round(P::Function, G::Function, M::Int, N::Int)
    cchain = CircuitGate{<:Any, M, <:Any}[]
    tmax = floor(Int, log(2, 2N/3))
    for t in tmax:-1:1
        mmax = floor(Int, (N-2^(t-1))/2^t)
        for m in 1:mmax
            push!(cchain, controlled_circuit_gate(( G(2^t*m) , P(t-1, 2m)), G(2^t*m + 2^(t-1)), X, M))
        end
    end
    return cchain
end

function c_round(P::Function, G::Function, M::Int, N::Int, Ñ::Int)
    cchain = CircuitGate{<:Any, M, <:Any}[]
    tmax = floor(Int, log(2, 2Ñ/3))
    for t in tmax:-1:1
        mmax = floor(Int, (Ñ-2^(t-1))/2^t)
        for m in 1:mmax
            push!(cchain, controlled_circuit_gate(( G(2^t*m) , P(t-1, 2m)), G(2^t*m + 2^(t-1)), X, M))
        end
    end
    return cchain
end

"""
    qcla_out_adder_circuit(N)

Specify wire functions for an out-of-place adder for 2 integers represented by `N` qubits.
Based on quantum carry-lookahead adder circuit on Draper et. al in quant-ph/0406142
"""

function qcla_out_adder_circuit(N)

    n = 1
    anc = 0
    while 2^n < N
        anc += N ÷ 2^n -1
        n += 1
    end
    M = 3N + anc + 1

    function a(m::Int)
        M - m
    end

    function b(m::Int)
        M - N - m
    end

    function s(m::Int)
        M - 2N - m
    end

    function G(m::Int)
        s(m)
    end

    function P(l::Int, m::Int)
        if l == 0
            return b(m)
        end

        while l > 1
            m += floor(Int, N/2^(l-1)-1)
            l -= 1
        end
        m
    end

    setup = CircuitGate{<:Any, M, <:Any}[]
    push!(setup, controlled_circuit_gate((a(0), b(0)), s(1), X, M))
    for i in 1:N-1
        push!(setup, controlled_circuit_gate((a(i), b(i)), s(i+1), X, M))
        push!(setup, controlled_circuit_gate(a(i), b(i), X, M))
    end
    cgc = CircuitGateChain{M}(setup)

    # Add gates for P rounds
    pchain = p_round(P, M, N)
    cgc *= CircuitGateChain{M}(pchain)

    # Add gates for G rounds
    gchain = g_round(P, G, M, N)
    cgc *= CircuitGateChain{M}(gchain)

    #  Add gates for C rounds
    cchain = c_round(P, G, M, N)
    cgc *= CircuitGateChain{M}(cchain)

    # Remove effect of P rounds
    cgc *= CircuitGateChain{M}(reverse(pchain))

    # Teardown
    teardown = CircuitGate{<:Any, M, <:Any}[]

    push!(teardown, controlled_circuit_gate(b(0), s(0), X, M))
    push!(teardown, controlled_circuit_gate(a(0), s(0), X, M))
    for i in 1:N-1
        push!(teardown, controlled_circuit_gate(b(i), s(i), X, M))
        push!(teardown, controlled_circuit_gate(a(i), b(i), X, M))
    end

    cgc *= CircuitGateChain{M}(teardown)
    return cgc
end


"""
    qcla_inplace_adder_circuit(N)

Specify wire functions for an out-of-place adder for 2 integers represented by `N` qubits.
Based on quantum carry-lookahead adder circuit on Draper et. al in quant-ph/0406142
"""

function qcla_inplace_adder_circuit(N)

    n = 1
    anc = 0
    while 2^n < N
        anc += N ÷ 2^n -1
        n += 1
    end
    M = 3N + anc

    function a(m::Int)
        m < N || error("a only takes m < N")
        M - m
    end

    function b(m::Int)
        m < N || error("b only takes m < N")
        M - N - m
    end

    function G(m::Int)
        m > 0 || error("G only takes positive integers")
        anc + m
    end

    function s(m::Int)
        m <= N || error("a only takes m <= N")
        if m < N
            return b(m)
        end
        G(m)
    end

    function P(l::Int, m::Int)
        if l == 0
            return b(m)
        end

        while l > 1
            m += floor(Int, N/2^(l-1)-1)
            l -= 1
        end
        m
    end

    # setup
    setup = CircuitGate{<:Any, M, <:Any}[]
    for i in 0:N-1
        push!(setup, controlled_circuit_gate((a(i), b(i)), G(i+1), X, M))
        push!(setup, controlled_circuit_gate(a(i), b(i), X, M))
    end
    setup = CircuitGateChain{M}(setup)

    # Add gates for P rounds
    pchain = p_round(P, M, N)
    cgc = CircuitGateChain{M}(pchain)

    # Add gates for G rounds
    gchain = g_round(P, G, M, N)
    cgc *= CircuitGateChain{M}(gchain)

    #  Add gates for C rounds
    cchain = c_round(P, G, M, N)
    cgc *= CircuitGateChain{M}(cchain)

    # Remove effect of P rounds
    cgc *= CircuitGateChain{M}(reverse(pchain))

    # Run intermediate circuit
    intermediate = CircuitGate{<:Any, M, <:Any}[]

    push!(intermediate, single_qubit_circuit_gate(b(0), X, M))
    for i in 1:N-2
        push!(intermediate, controlled_circuit_gate(G(i), b(i), X, M))
        push!(intermediate, single_qubit_circuit_gate(b(i), X, M))
        push!(intermediate, controlled_circuit_gate(a(i), b(i), X, M))
    end
    if N > 1
        push!(intermediate, controlled_circuit_gate(G(N-1), b(N-1), X, M))
    end

    if N > 1
        pchain2 = p_round(P, M, N, N-1)
        gchain2 = g_round(P, G, M, N, N-1)
        cchain2 = c_round(P, G, M, N, N-1)
        rpchain2 = reverse(pchain2)
        r = vcat(pchain2, gchain2, cchain2, rpchain2)
    else
        r = CircuitGate{<:Any, M, <:Any}[]
    end
    cgc *= CircuitGateChain{M}(intermediate) * CircuitGateChain{M}(reverse(r))

    # ## Teardown
    teardown = CircuitGate{<:Any, M, <:Any}[]

    push!(teardown, controlled_circuit_gate( (a(0), b(0)), G(1), X, M))
    push!(teardown, single_qubit_circuit_gate(b(0), X, M))
    for i in 1:N-2
        push!(teardown, controlled_circuit_gate(a(i), b(i), X, M))
        push!(teardown, controlled_circuit_gate( (a(i), b(i)), G(i+1), X, M))
        push!(teardown, single_qubit_circuit_gate(b(i), X, M))
    end


    cgc = setup * cgc * CircuitGateChain{M}(teardown)
    return cgc
end

# """
#     qcla_comparator_circuit(N)
#
# Create circuit for QCLA comparator for 2 integers represented by `N` qubits.
# Based on quantum carry-lookahead comparator circuit on Draper et. al in quant-ph/0406142
# """
#
# function qcla_comparator_circuit(N)
#
#     n = 1
#     anc = N - floor(Int, log(2, N-1)) - 2
#     M = 3N + anc
#
#     function a(m::Int)
#         m < N || error("a only takes m < N")
#         M - m
#     end
#
#     function b(m::Int)
#         m < N || error("b only takes m < N")
#         M - N - m
#     end
#
#     function G(m::Int)
#         m > 0 || error("G only takes positive integers")
#         anc + m
#     end
#
#     function s(m::Int)
#         m <= N || error("a only takes m <= N")
#         if m < N
#             return b(m)
#         end
#         G(m)
#     end
#
#     function P(l::Int, m::Int)
#         if l == 0
#             return b(m)
#         end
#
#         while l > 1
#             m += floor(Int, N/2^(l-1)-1)
#             l -= 1
#         end
#         m
#     end
#
#     # setup
#     setup = CircuitGate{<:Any, M, <:Any}[]
#         push!(setup, single_qubit_circuit_gate(a(0), X, M))
#         push!(setup, controlled_circuit_gate((a(0), b(0)), G(1), X, M))
#     for i in 1:N-1
#         push!(setup, single_qubit_circuit_gate(a(i), X, M))
#         push!(setup, controlled_circuit_gate((a(i), b(i)), G(i+1), X, M))
#         push!(setup, controlled_circuit_gate(a(i), b(i), X, M))
#     end
#
#     setup = CircuitGateChain{M}(setup)
#
#     # Add gates for P rounds
#     pchain = p_round(P, M, N)
#     cgc = CircuitGateChain{M}(pchain)
#
#     # Add gates for G rounds
#     gchain = g_round(P, G, M, N)
#     cgc *= CircuitGateChain{M}(gchain)
#
#     if N > 1
#         gchain2 = g_round(P, G, M, N, N-1)
#     else
#         gchain2 = CircuitGate{<:Any, M, <:Any}[]
#     end
#
#     cgc *= CircuitGateChain{M}(reverse(gchain2)) * CircuitGateChain{M}(reverse(pchain))
#
#     # Run intermediate circuit
#     teardown = CircuitGate{<:Any, M, <:Any}[]
#     push!(teardown, controlled_circuit_gate((a(0), b(0)), G(1), X, M))
#     push!(teardown, single_qubit_circuit_gate(a(0), X, M))
#     for i in 1:N-2
#         push!(teardown, controlled_circuit_gate((a(i), b(i)), G(i+1), X, M))
#         push!(teardown, controlled_circuit_gate(a(i), b(i), X, M))
#         push!(teardown, single_qubit_circuit_gate(a(i), X, M))
#     end
#     if N > 1
#         push!(teardown, controlled_circuit_gate(a(N-1), b(N-1), X, M))
#         push!(teardown, single_qubit_circuit_gate(a(N-1), X, M))
#     end
#     push!(teardown, single_qubit_circuit_gate(G(N), X, M))
#
#     cgc = setup * cgc * CircuitGateChain{M}(teardown)
#     return cgc
# end
