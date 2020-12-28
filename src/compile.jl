using Memoize
using GenericSchur
using JLD

"""
    greyencode(n::Integer) = n ⊻ (n >> 1)

returns grey encoding of given integer
"""
greyencode(n::Integer) = n ⊻ (n >> 1)

"""
    inverseM(N)

finds inverse of M matrix of size 2^(N-1)×2^(N-1) with elements M_{i,j} equal
to inner product of binary rep of i (classical) and j (grey code), i,j ∈ [1:2^(N-1)]
utilizes algorithm from arxiv:quant-ph/0404089
"""
@memoize function inverseM(N)
    m = 2^(N - 1)
    M = zeros((m, m))
    for i in 1:m
        for j in 1:m
            M[i,j] = (-1)^count_ones((i - 1) .& greyencode(j - 1))
        end
    end
    return inv(M)
end

"""
    stateprep(u, N)

creates Array of CircuitGate objects preparing state `u` over `N` qubits..
uses algorithm from arxiv:quant-ph/0410066
"""
function stateprep(u, N, n=1)
    cg = CircuitGate[]
    k = 2
    grey = zeros(Int, 2^(N - n))
    for i in 0:N - 2
        grey[2^i:2^(i + 1):end] .= i + n + 1
    end
    grey[end] = N
    θ = real.(atan.(-u[2:2:end] ./ u[1:2:end]) .* 2)
    θ[isnan.(θ)] .= -π
    M = inverseM(N - n + 1)
    θ = M * θ
    for i in 1:2^(N - n)
        push!(cg, circuit_gate((n,), RyGate(θ[i])))
        push!(cg, circuit_gate((n), X, grey[i]))
    end
    return cg
end

"""
    compile!(m::Matrix{Complex64,2}) where {N}

compiles quantum circuit from given unitary matrix. usage of algorithms from
ff10.1016/j.cpc.2019.107001f, 10.1103/PhysRevA.69.032315, arxiv.org:1003.5760, arxiv:quant-ph/0404089
"""
function compile(m::AbstractMatrix{ComplexF64}, N, wires=nothing)
    m * m' ≈ I || error("Only unitary matrices can be compiled into a quantum circuit")
    s = size(m)
    s[1] == s[2] || error("Only square matrices can be compiled into a quantum circuit")
    
    if Diagonal(diag(m)) ≈ m
        return compilediag(diag(m), N, wires)
    end
    
    if size(m)[1] == 2
        return compile1qubit(m, wires)[1]
    end
    
    if size(m)[1] == 4
        return compile2qubit(m, N, wires)
    end
    
    cgc = CircuitGate[]
    QR, τ = qr_unblocked(m)
    cgc = compilediag(diag(QR), N, wires)
    
    Dg = fill(1.0 + 0.0im, (2^N))
    Dg[1] = -1
    Dgm = compilediag(Dg, N)
    
    for i in 2^(N) - 1:-1:1
        prepd = CircuitGate[]
        u = zeros(ComplexF64, (2^N))
        u[i + 1:end] = QR[i + 1:end, i]
        u[i] = 1
        u = u ./ norm(u)
        angles = exp.(im .* angle.(u))
        u = u ./ angles
        D = compilediag(angles, N)

        for j in 1:N - 1
            cg = stateprep(u[1:2^(j - 1):end], N, j)
            append!(prepd, cg)
            u = apply(cg, u)
        end
        θ = real(atan(-u[2^(N - 1) + 1] ./ u[1]) .* 2)
        if !isnan(θ)
            cg = circuit_gate((N,), RyGate(θ))
        else
            cg = circuit_gate((N,), RyGate(-π))
        end
        push!(prepd, cg)

        prep = reverse(adjoint.(prepd))
        append!(prep, D)
        prepend!(prepd, adjoint.(D))

        append!(cgc, prepd)
        append!(cgc, Dgm)
        append!(cgc, prep)
    end
    return cgc
end

"""
    qr_blocked(m::Matrix{Complex64,2})

performs blocked qr decomposition on matrix `m`
"""
function qr_blocked!(m::AbstractMatrix{ComplexF64}, block_size=2::Integer)
    N = size(m)[1]
    N == size(m)[2] || error("Matrix `m` must be a square matrix")
    if block_size > N
        return qr_unblocked(m)
    end
    τ = ComplexF64[]
    num_blocks = N ÷ block_size
    for x in 0:num_blocks - 1
        row = x * block_size
        col = x * block_size
        for y in 1:block_size
                a = angle(m[row + y, col + y])
                b = norm(m[row + y, col + y])
                m[row + y,col + y] = -exp(im * a)
                m[row + y + 1:N, col + y] =  m[row + y + 1:N, col + y] ./ (exp(im * a) * (1 + b))
                push!(τ, 1 + b)
            if y < block_size
                update1 = m[row + y + 1:N, col + y] .* m[row + y, col + y + 1]
                update2 = m[row + y, col + y + 1:N] .* m[row + y + 1, col + y]
                m[row + y + 1:N,col + y + 1] = m[row + y + 1:N,col + y + 1] - update1
                m[row + y + 1,col + y + 2:N] = m[row + y + 1,col + y + 2:N] - update2[2:end]
            end
        end
        if col + block_size < N
            m[row + block_size + 1:N, col + block_size + 1:N] -=
                m[row + block_size + 1:N, col + 1:col + block_size] *
                m[row + 1:row + block_size, col + block_size + 1:N]
        end
    end
    if N % block_size != 0
        m_unblocked, τ_unblocked = qr_unblocked(m[num_blocks * block_size:end, num_blocks * block_size:end])
        m[num_blocks * block_size:end,num_blocks * block_size:end] = m_unblocked
        append!(τ, τ_unblocked)
    end
    return m, τ
end

"""
    qr_unblocked(m::Matrix{Complex64,2})

performs serial qr decomposition on matrix `m`
"""
function qr_unblocked!(m::AbstractMatrix{ComplexF64}, n=1::Integer)
    N = size(m)[1]
    τ = ComplexF64[]
    for i in n:N - 1
        b = norm(m[i,i])
        a = angle(m[i,i])
        push!(τ, 1 + b)
        m[i,i] = -exp(im * a)
        m[i + 1:N, i] =  m[i + 1:N, i] ./ (exp(im * a) * (1 + b))
        m[i + 1:N,i + 1:N] -= (m[i + 1:N, i] * transpose(m[i,i + 1:N]))
    end
    return τ
end

"""
    qr_unblocked(m::Matrix{Complex64,2})

performs serial qr decomposition on matrix `m`
"""
function qr_unblocked(m::AbstractMatrix{ComplexF64}, n=1::Integer)
    m = deepcopy(m)
    N = size(m)[1]
    τ = ComplexF64[]
    for i in n:N - 1
        b = norm(m[i,i])
        a = angle(m[i,i])
        push!(τ, 1 + b)
        m[i,i] = -exp(im * a)
        m[i + 1:N, i] =  m[i + 1:N, i] ./ (exp(im * a) * (1 + b))
        m[i + 1:N,i + 1:N] -= (m[i + 1:N, i] * transpose(m[i,i + 1:N]))
    end
    return m, τ
end

"""
    compile1qubit(m::AbstractMatrix{ComplexF64}, N, wires=nothing)

compiles a U(2) matrix into quantum gates
"""
function compile1qubit(m::AbstractMatrix{ComplexF64}, wires=nothing)
    phase = sqrt(m[1,1] * m[2,2] / norm(m[1,1] * m[2,2]))
    b = m / phase
    ϕ = real(acos(sqrt(b[1,1] * b[2,2])) * 2)
    θ1 = imag(log(2im * b[2,1] * b[2,2] / sin(ϕ)))
    
    θ2 = imag(-log(2 * im * b[1,1] * b[2,1] / sin(ϕ)))
    if isnothing(wires)
        wires = [1]
    end
    cg = CircuitGate[]
    append!(cg, circuit_gate.((wires[1],), [RzGate(θ2), RxGate(ϕ), RzGate(θ1)]))
    return cg, phase
end

"""
    decomposeSO4(m::AbstractMatrix{ComplexF64})

decomposes a SO(4) matrix into 2 SU(2) matrices, such that `m`∈SO(4), A,B ∈ SU(2)
via isoclinic decomposition. m ≊ A ⊗ B under the magic basis homomorphism.
"""
function decomposeSO4(m::AbstractMatrix{ComplexF64})
    size(m)[1] == size(m)[2] || error("decomposeSO4 only works on square matrices")
    size(m)[1] == 4 || error("decomposeSO4 only works on 4x4 matrices")

    ap = (m[1,1] + m[2,2] + m[3,3] + m[4,4]) / 4
    bp = (m[2,1] - m[1,2] + m[4,3] - m[3,4]) / 4
    cp = (m[3,1] - m[4,2] - m[1,3] + m[2,4]) / 4
    dp = (m[4,1] + m[3,2] - m[2,3] - m[1,4]) / 4

    aq = (m[2,1] - m[1,2] - m[4,3] + m[3,4]) / 4
    ar = (m[3,1] + m[4,2] - m[1,3] - m[2,4]) / 4
    as = (m[4,1] - m[3,2] + m[2,3] - m[1,4]) / 4

    p = real(sqrt(ap^2 + bp^2 + cp^2 + dp^2))

    a = real(ap / p)
    b = real(bp / p)
    c = real(cp / p)
    d = real(dp / p)

    q = real(aq / a)
    r = real(ar / a)
    s = real(as / a)

    if (p^2 + q^2 + r^2 + s^2) ≈ 1
        A = [p + q * im s + r * im; -s + r * im p - q * im]
        B = [a + b * im -d + c * im; d + c * im a - b * im]
        return A, B
    end

    p = -p

    a = real(ap / p)
    b = real(bp / p)
    c = real(cp / p)
    d = real(dp / p)

    q = real(aq / a)
    r = real(ar / a)
    s = real(as / a)

    @assert  (p^2 + q^2 + r^2 + s^2) ≈ 1

    A = [p + q * im s + r * im; -s + r * im p - q * im]
    B = [a + b * im -d + c * im; d + c * im a - b * im]

    return A, B
end

"""
    compile2qubit(m::AbstractMatrix{ComplexF64}, N, wires=nothing)

compiles an arbitrary U(4) matrix into a quantum circuit. Algorithm taken from
arxiv:quant-ph/0308006, arxiv:quant-ph/0211002
"""
function compile2qubit(m::AbstractMatrix{ComplexF64}, N, wires=nothing)
    cg = CircuitGate[]
    E = 1 / sqrt(2) .* [1 im 0 0; 0 0 im 1; 0 0 im -1; 1 -im 0 0]

    U = E' * m * E
    P2 = U * transpose(U)

    Diag, K_2 = eigen(P2)
    Diag[1] = Diag[1] / det(K_2)^2
    K_2[:,1] = K_2[:,1] * det(K_2)

    Diag = diagm(sqrt.(Diag))
    Diag[1] = Diag[1] * det(U) / det(Diag)

    P = K_2 * Diag * inv(K_2)
    K_1 = inv(P) * U
    C, D = decomposeSO4(inv(K_2) * K_1)

    append!(cg, compile1qubit(C, [2])[1])
    append!(cg, compile1qubit(D, [1])[1])

    append!(cg, [circuit_gate(2, X, 1),
                    circuit_gate((1,), HadamardGate()),
                    circuit_gate((2,), SdagGate()),
                    circuit_gate((1,), SdagGate())])
    append!(cg, compilediag(diag(Diag), N))
    append!(cg, [circuit_gate((1,), SGate()),
                    circuit_gate((2,), SGate()),
                    circuit_gate((1,), HadamardGate()),
                    circuit_gate((2), X, 1)])

    A, B = decomposeSO4(K_2)

    append!(cg, compile1qubit(A, [2])[1])
    append!(cg, compile1qubit(B, [1])[1])

    return cg
end

"""
    fillψ!(d::Vector{ComplexF64}, ψ::Vector{Float64}, l::Integer)

fills empty vector `ψ`, such that  ψ = −im*[logχ1(d) log χ2(d) ··· logχl−1(U)],
where χn(d) = d[2n-1]*d[2n+2]/(d[2n]*d[2n+1]) per arxiv:quant-ph/0303039
"""
function fillψ!(d::Vector{ComplexF64}, ψ::Vector{Float64}, l::Integer)
    for i in StepRange(1, 1, l - 1)
        ψ[i] = imag(log(d[2i - 1] * d[2i + 2] / (d[2i] * d[2i + 1])))
    end
end

"""
    ηcol(l::Integer, n::Integer)

calculates a single flip state for given integer `n` and a column length of `l`,
    where 0 <= n <= l. returns a vector of length l, where vec[n] = 1, vec[n+1] = -1.
    if n == 0, vec[1] = -1, if n==l, vec[l] = 1
# Examples
```julia-repl
julia> l = 4
julia> ηcol(4, 0)
[-1, 0, 0, 0]
julia> ηcol(4, 1)
[1, -1, 0, 0]
julia> ηcol(4, 2)
[0, 1, -1, 0]
julia> ηcol(4, 3)
[0, 0, 1, -1]
julia> ηcol(4, 4)
[0, 0, 0, 1]
```
"""
@memoize Dict function ηcol(l::Integer, n::Integer)
    vec = zeros(Integer, l)
    if n == 0
        vec[1] = -1
        return vec
    elseif n == l
    vec[l] = 1
        return vec
    end

    vec[n] = 1
    vec[n + 1] = -1
    vec
end

"""
    flip_state(m::Integer, l::Integer)

calculates the `m`th column of the η⊗ matrix by summing all flip states
associated with the binary representation of `m` over the span of [1,l-1]
algorithm taken from arxiv:quant-ph/0303039
"""
@memoize Dict function flip_state(m::Integer, l::Integer)
    # calculates flip states of a given integer `x`
    f = filter(x -> count_ones(x & m) % 2 == 1, 1:l - 1)
    sum(ηcol.((l - 1,), f))
end

"""
    svalue(i::Integer)

returns tuple of position of 1s in the binary representation of integer `i` starting from least significant digit.

# Examples
```julia-repl
julia> svalue(2)
(2,)
julia> svalue(4)
(3,)
julia> svalue(7)
(1,2,3)
```
"""
function svalue(i::Integer)
    b = Int[]
    count = 1
    while i != 0
        if i & 1 == 1
            push!(b, count)
        end
        count += 1
i = i >> 1
    end
    Tuple(b)
end

"""
    compilediag(d::Vector{ComplexF64}, N, cg=nothing, j=0)

recursively compiles an arbitrary diagonal unitary matrix D_{2^N} ∈ C^(2^N×2^N) into a quantum circuit. algorithm taken from arxiv:quant-ph/0303039
D_{2^N} can be decomposed into D_{2^(N-1)} ∈ C^(2^(N-1)×2^(N-1)) ⊗ exp(i*ϕ)*RzGate(θ)
"""
function compilediag(d::Vector{ComplexF64}, N, cg=nothing, j=0)
    cgs = CircuitGate[]
    # single qubit diagonal matrix can be split into a phase shift + RzGate
    if length(d) == 2
        logd = imag.(log.(d))
        β = (logd[1] - logd[2])
        if isnothing(cg)
            cg = CircuitGate[]
        end
        push!(cg, circuit_gate((j + 1,), RzGate(-β)))
        return cg
    end
    l = 2^(N - j - 1)
    ψ = zeros(Float64, l - 1)
    α = zeros(Float64, l)

    fillψ!(d, ψ, l)

    ηplus = zeros(Float64, (l - 1, l - 1))
    for i in StepRange(1, 1, l - 1)
    ηplus[:,i] = flip_state(greyencode(i), l)
    end

    ηplus = inv(ηplus)
    α[2:end] = -0.5 * ηplus * ψ

    gatewires = svalue(greyencode(1))
    oldwires = ()

    for i in StepRange(2, 1, l)
        rz = [exp(-im * α[i] / 2), exp(im * α[i] / 2)]
        for _ in 1:N - j - 1
            rz = append!(rz, rz)
        end
        s = filter(x -> count_ones(x & greyencode(i - 1)) % 2 == 1, 1:l - 1)
        s = append!(2 .* s .+ 1, 2 .* s .+ 2)

        for k in symdiff(gatewires, oldwires)
            push!(cgs, circuit_gate((j + 1), X, k + 1 + j))
        end
        oldwires = gatewires
        gatewires = svalue(greyencode(i))
        push!(cgs, circuit_gate((j + 1,), RzGate(α[i])))

        rz[s] .= 1 ./ rz[s]
        d = d ./ rz
    end

    for k in reverse(oldwires)
        push!(cgs, circuit_gate(j + 1, X, k + 1 + j))
    end

    if !(d[1] ≈ d[2])
        logd = imag.(log.(d))
        α[1] = logd[1] - logd[2]
        pushfirst!(cgs, circuit_gate((j + 1,), RzGate(-α[1])))
        d[1:2:end] = d[1:2:end] ./ exp(im * α[1] / 2)
    end
    if isnothing(cg)
        return compilediag(d[1:2:end], N, cgs, j + 1)
    end
    return compilediag(d[1:2:end], N, append!(cg, cgs), j + 1)
end
