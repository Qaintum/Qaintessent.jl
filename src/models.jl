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
    toffoli_circuit(cntrl, trg, N)

construct the circuit for decomposing the Toffoli gate in Figure 4.9 of Nielsen and Chuang (2000). the constructed toffoli gate acts in circuit of `N` qubits, has controls
    on wires in Tuple `cntrl` and has target on wire `trg`.
returns a `CircuitGateChain{N}` object.
"""
function toffoli_circuit(cntrl::Tuple{<:Integer, <:Integer} , trg::Integer, N::Integer)
    CircuitGateChain{N}([
        single_qubit_circuit_gate(trg,                HadamardGate(), N),
        controlled_circuit_gate(  cntrl[2], trg,      X,              N),
        single_qubit_circuit_gate(trg,                TdagGate(),     N),
        controlled_circuit_gate(  cntrl[1], trg,      X,              N),
        single_qubit_circuit_gate(trg,                TGate(),        N),
        controlled_circuit_gate(  cntrl[2], trg,      X,              N),
        single_qubit_circuit_gate(trg,                TdagGate(),     N),
        controlled_circuit_gate(  cntrl[1], trg,      X,              N),
        single_qubit_circuit_gate(trg,                TGate(),        N),
        single_qubit_circuit_gate(cntrl[2],           TGate(),        N),
        controlled_circuit_gate(  cntrl[1], cntrl[2], X,              N),
        single_qubit_circuit_gate(trg,                HadamardGate(), N),
        single_qubit_circuit_gate(cntrl[1],           TGate(),        N),
        single_qubit_circuit_gate(cntrl[2],           TdagGate(),     N),
        controlled_circuit_gate(  cntrl[1], cntrl[2], X,              N),
    ])
end


```jldoctest
N = 4
adder = vbe_adder_circuit(N)
ψ = fill(0.0+0.0*im, 2^M)
a = 3
b = 5

index = b << N + a
ψ[index+1] = 1.0

ψ = apply(adder, ψ)
((findall(x->x==1, ψ)[1]-1)%(2^2N)) >> N

# output

8
```

"""
    vbe_adder_circuit(N)

Construct an in-place adder for 2 integers represented by `N` qubits.
Based on ripple-carry adder circuit in Vedral et. al (Phys. Rev. A 54, 147 (1996), arXiv:quant-ph/9511018)
Returns a CircuitGateChain{3N+1} as there are N+1 ancillary wires

If the two added integers are represented as:

    ``A = a_{0}\\times 2^{0} + a_{1} \\times 2^{1} + a_{2} \\times 2^{2} + .. + a_{N} \\times 2^{N} \\\\ B = b_{0}\\times 2^{0} + b_{1} \\times 2^{1} + b_{2} \\times 2^{2} + .. + b_{N} \\times 2^{N}``

The input index should be ``a_{0}a_{1}a_{2}..a_{N}b_{0}b_{1}b_{2}..b_{N} + 1`` with ``a_{0}`` as the fastest running index

i.e.

    input = fill(0, 2^(3N+1))
    index = b << N + a
    input[index+1] = 1

The output index will be in the form ``a_{0}a_{1}a_{2}..a_{N}c_{0}c_{1}c_{2}..c_{N} + 1`` where:

    ``C = (A+B) % (2^{N} + 1) = c_{0}\\times 2^{0} + c_{1} \\times 2^{1} + c_{2} \\times 2^{2} + ... + c_{N} \\times 2^{N}``

i.e.

    output = apply(adder, input)
    answer = ((findall(x->x==1, output)[1]-1)%(2^2N)) >> N
    @ test answer == (a+b)%(2^N)
"""
function vbe_adder_circuit(N::Integer)
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
function p_round(P::Function, M::Integer, N::Integer)
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

function p_round(P::Function, M::Integer, N::Integer, Ñ::Integer)
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
function g_round(P::Function, G::Function, M::Integer, N::Integer)
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

function g_round(P::Function, G::Function, M::Integer, N::Integer, Ñ::Integer)
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
function c_round(P::Function, G::Function, M::Integer, N::Integer)
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

function c_round(P::Function, G::Function, M::Integer, N::Integer, Ñ::Integer)
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

```jldoctest
N = 4
adder = qcla_out_adder_circuit(N)
ψ = fill(0.0+0.0*im, 2^M)
a = 2
b = 6

index = b << N + a
ψ = fill(0.0+0.0*im, 2^M)

ψ[index+1] = 1.0

ψ = apply(cgc, ψ)
(findall(x->x==1, ψ)[1] - 1) >> 2N

# output

8
```

"""
    qcla_out_adder_circuit(N)

Construct an out-of-place adder for 2 integers represented by `N` qubits. returns a `CircuitGateChain{3N+1}` object.
Based on quantum carry-lookahead adder circuit by Draper et. al (Quant. Inf. Comp. 6, 351-369 (2006), arXiv:quant-ph/0406142)
Returns a CircuitGateChain{3N+1} as there are N+1 ancillary wires

If the two added integers are represented as:

    ``A = a_{0}\\times 2^{0} + a_{1} \\times 2^{1} + a_{2} \\times 2^{2} + .. + a_{N} \\times 2^{N} \\\\ B = b_{0}\\times 2^{0} + b_{1} \\times 2^{1} + b_{2} \\times 2^{2} + .. + b_{N} \\times 2^{N}``

The input index should be ``a_{0}a_{1}a_{2}..a_{N}b_{0}b_{1}b_{2}..b_{N} + 1`` as Julia starts indexing at `1` and ``a_{0}`` as the fastest running index

i.e.

    input = fill(0, 2^(3N+1))
    index = b << N + a
    input[index+1] = 1

The output index will be in the form ``a_{0}a_{1}a_{2}..a_{N}b_{0}b_{1}b_{2}..b_{N}c_{0}c_{1}c_{2}..c_{N+1} + 1`` where:

    ``C = A+B = c_{0}\\times 2^{0} + c_{1} \\times 2^{1} + c_{2} \\times 2^{2} + ... + c_{N+1} \\times 2^{N+1}``

i.e.

    output = apply(adder, input)
    answer = (findall(x->x==1, ψ)[1] - 1) >> 2N
    @ test answer == a+b
"""
function qcla_out_adder_circuit(N)

    n = 1
    anc = 0
    while 2^n < N
        anc += N ÷ 2^n -1
        n += 1
    end
    M = 3N + anc + 1

    function a(m::Integer)
        M - m
    end

    function b(m::Integer)
        M - N - m
    end

    function s(m::Integer)
        M - 2N - m
    end

    function G(m::Integer)
        s(m)
    end

    function P(l::Integer, m::Integer)
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

```jldoctest
N = 4
adder = qcla_inplace_adder_circuit(N)
ψ = fill(0.0+0.0*im, 2^M)
a = 6
b = 4

index = b << N + a
ψ = fill(0.0+0.0*im, 2^M)

ψ[index+1] = 1.0

ψ = apply(cgc, ψ)
(findall(x->x==1, ψ)[1] - 1) >> N

# output

10
```


"""
    qcla_inplace_adder_circuit(N)

Construct an in-place adder for 2 integers represented by `N` qubits.
Based on quantum carry-lookahead adder circuit by Draper et. al (Quant. Inf. Comp. 6, 351-369 (2006), arXiv:quant-ph/0406142)
Returns a CircuitGateChain{3N+1} as there are N+1 ancillary wires

If the two added integers are represented as:

    ``A = a_{0}\\times 2^{0} + a_{1} \\times 2^{1} + a_{2} \\times 2^{2} + .. + a_{N} \\times 2^{N} \\\\ B = b_{0}\\times 2^{0} + b_{1} \\times 2^{1} + b_{2} \\times 2^{2} + .. + b_{N} \\times 2^{N}``

The input index should be ``a_{0}a_{1}a_{2}..a_{N}b_{0}b_{1}b_{2}..b_{N} + 1`` as Julia starts indexing at `1` and ``a_{0}`` as the fastest running index

i.e.

    input = fill(0, 2^(3N+1))
    index = b << N + a
    input[index+1] = 1

The output index will be in the form ``a_{0}a_{1}a_{2}..a_{N}b_{0}c_{1}c_{2}..c_{N+1} + 1`` where:

    ``C = A+B = c_{0}\\times 2^{0} + c_{1} \\times 2^{1} + c_{2} \\times 2^{2} + ... + c_{N+1} \\times 2^{N+1}``

i.e.

    output = apply(adder, input)
    answer = (findall(x->x==1, ψ)[1] - 1) >> N
    @ test answer == a+b
"""
function qcla_inplace_adder_circuit(N)

    n = 1
    anc = 0
    while 2^n < N
        anc += N ÷ 2^n -1
        n += 1
    end
    M = 3N + anc

    function a(m::Integer)
        m < N || error("a only takes m < N")
        M - m
    end

    function b(m::Integer)
        m < N || error("b only takes m < N")
        M - N - m
    end

    function G(m::Integer)
        m > 0 || error("G only takes positive integers")
        anc + m
    end

    function s(m::Integer)
        m <= N || error("a only takes m <= N")
        if m < N
            return b(m)
        end
        G(m)
    end

    function P(l::Integer, m::Integer)
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

    if N > 1
        push!(teardown, controlled_circuit_gate( (a(0), b(0)), G(1), X, M))
    end
    push!(teardown, single_qubit_circuit_gate(b(0), X, M))
    for i in 1:N-2
        push!(teardown, controlled_circuit_gate(a(i), b(i), X, M))
        push!(teardown, controlled_circuit_gate( (a(i), b(i)), G(i+1), X, M))
        push!(teardown, single_qubit_circuit_gate(b(i), X, M))
    end


    cgc = setup * cgc * CircuitGateChain{M}(teardown)
    return cgc
end

# TODO: Add comparator circuit for QCLA adder
