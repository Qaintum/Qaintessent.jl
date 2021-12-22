using Memoize
using LinearAlgebra

"""
    greyencode(n::Int) = n ⊻ (n >> 1)
returns grey encoding of given Int
"""
greyencode(n::Int) = n ⊻ (n >> 1)

"""
    inverseM(N)
finds inverse of M matrix of size 2^(N-1)×2^(N-1) with elements M_{i,j} equal to
inner product of binary rep of i (classical) and j (grey code), i,j ∈
[1:2^(N-1)] utilizes algorithm from arxiv:quant-ph/0404089
"""
@memoize function inverseM(N::Int)
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
    
creates Array of CircuitGate objects preparing state `u` over `N` qubits.. uses
algorithm from arxiv:quant-ph/0410066
"""
function stateprep(u::Vector{<:Complex}, N::Int, n::Int=1)
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
    unitary2circuit(m::Matrix{Complex64,2}) where {N}

compiles quantum circuit from given unitary matrix. Implements algorithms from
ff10.1016/j.cpc.2019.107001f, 10.1103/PhysRevA.69.032315, arxiv.org:1003.5760,
arxiv:quant-ph/0404089

# Examples
```@example; setup = :(using RandomMatrices)
M = Stewart(ComplexF64, 4);
cgc = unitary2circuit(M)
```
"""
function unitary2circuit(m::AbstractMatrix{ComplexF32}, N::Union{Int, Nothing}=nothing, wires=nothing)
    unitary2circuit(convert(Matrix{ComplexF64}, m), N, wires)
end

function unitary2circuit(m::AbstractMatrix{ComplexF64}, N::Union{Int, Nothing}=nothing, wires=nothing)
    isapprox(m * m', I, atol=1e2*eps(FloatQ)) || error("Only unitary matrices can be compiled into a quantum circuit")
    s = size(m)
    
    if isnothing(N)
        N = intlog2(s[1])
    end
    
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
        u = zeros(ComplexQ, (2^N))
        u[i + 1:end] = QR[i + 1:end, i]
        u[i] = 1
        u = u ./ norm(u)
        angles = exp.(im .* angle.(u))
        u = u ./ angles
        D = compilediag(angles, N)

        for j in 1:N - 1
            cg = stateprep(u[1:2^(j - 1):end], N, j)
            append!(prepd, cg)
            u = apply(u, cg)
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
    qr_unblocked(m::Matrix{Complex64,2})

performs serial qr decomposition on matrix `m`
"""
function qr_unblocked!(m::AbstractMatrix{ComplexF64}, n=1::Int)
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
function qr_unblocked(m::AbstractMatrix{ComplexF64}, n=1::Int)
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
function compile1qubit(m::AbstractMatrix{<:Complex}, wires=nothing)
    if isnothing(wires)
        wires = [1]
    end
    cg = CircuitGate[]
    phase = sqrt(m[1,1]*m[2,2] - m[2,1]*m[1,2])
    m  = m / phase
    b = norm.(m)
    if (norm(b[1,1]) + norm(b[2,2])) < 1e-5
        θ2 = -imag(log(m[2,1]/m[1,2]))
        θ1 = 0
        ϕ = (real(m[1,2]/exp(-θ1/2*im)) < 0 ? 3π/2 : π/2) * 2
    elseif (norm(b[1,2]) + norm(b[2,1])) < 1e-5
        θ2 = -imag(log(m[1,1]/m[2,2]))
        θ1 = 0
        ϕ = (real(m[1,1]/exp(-θ1/2*im)) < 0 ? π : 0) * 2
    else
        ϕ = real(acos(sqrt(m[1,1] * m[2,2])) * 2)
        θ1 = imag(log(2im * m[2,1] * m[2,2] / sin(ϕ)))
        θ2 = imag(-log(2im * m[1,1] * m[2,1] / sin(ϕ)))
    end

    if !(norm(θ2) < 1e-5)
        append!(cg, circuit_gate.((wires[1],), [RzGate(θ2), RxGate(ϕ)]))
    else
        append!(cg, circuit_gate.((wires[1],), [RxGate(ϕ)]))
    end

    if !(norm(θ1) < 1e-5)
        append!(cg, circuit_gate.((wires[1],), [RzGate(θ1)]))
    end
    return cg, phase
end

iso_permutations = [
                    [[1,1], [2,2], [3,3], [4,4]],
                    [[2,1], [1,2], [4,3], [3,4]],
                    [[3,1], [4,2], [1,3], [2,4]],
                    [[4,1], [3,2], [2,3], [1,4]],
                    [[2,1], [1,2], [4,3], [3,4]],
                    [[1,1], [2,2], [3,3], [4,4]],
                    [[4,1], [3,2], [2,3], [1,4]],
                    [[3,1], [4,2], [1,3], [2,4]],
                    [[3,1], [4,2], [1,3], [2,4]],
                    [[4,1], [3,2], [2,3], [1,4]],
                    [[1,1], [2,2], [3,3], [4,4]],
                    [[2,1], [1,2], [4,3], [3,4]],
                    [[4,1], [3,2], [2,3], [1,4]],
                    [[3,1], [4,2], [1,3], [2,4]],
                    [[2,1], [1,2], [4,3], [3,4]],
                    [[1,1], [2,2], [3,3], [4,4]],
                    ]

iso_signs = [
            [ 1, 1, 1, 1],
            [ 1,-1, 1,-1],
            [ 1,-1,-1, 1],
            [ 1, 1,-1,-1],
            [ 1,-1,-1, 1],
            [-1,-1, 1, 1],
            [-1,-1,-1,-1],
            [ 1,-1, 1,-1],
            [ 1, 1,-1,-1],
            [ 1,-1,-1, 1],
            [-1, 1,-1, 1],
            [-1,-1,-1,-1],
            [ 1,-1, 1,-1],
            [-1,-1,-1,-1],
            [ 1, 1,-1,-1],
            [-1, 1, 1,-1]
            ]
            

function _get_permutation(m::AbstractMatrix{Float64}, ::Type{Val{I}}, ::Type{Val{J}}) where {I,J}
    val = 0    
    for i in 1:4
        val += m[iso_permutations[4*(J-1)+I][i]...] * iso_signs[4*(J-1)+I][i]
    end
    val/4
end

function _get_col!(m::AbstractMatrix{Float64}, col::Vector{Float64}, col_id::Int)
    for i in 1:4
        col[i] = _get_permutation(m, Val{i}, Val{col_id})
    end
end

function _get_row!(m::AbstractMatrix{Float64}, row::Vector{Float64}, row_id::Int)
    for i in 1:4
        row[i] = _get_permutation(m, Val{row_id}, Val{i})
    end
end

"""
    decomposeSO4(m::AbstractMatrix{ComplexF64})

decomposes a SO(4) matrix into 2 SU(2) matrices, such that `m`∈SO(4), A,B ∈
SU(2) via isoclinic decomposition. m ≊ A ⊗ B under the magic basis homomorphism.
"""
function decomposeSO4(m::AbstractMatrix{Float64})
    size(m)[1] == size(m)[2] || error("decomposeSO4 only works on square matrices")
    size(m)[1] == 4 || error("decomposeSO4 only works on 4x4 matrices")
    isapprox(det(m), 1; rtol=1e-5) || error("matrix `m` is not in SO(4), determinant is not 1")

    col = zeros(Float64, 4)
    row = zeros(Float64, 4)

    col_id = 0
    row_id = 0

    for i in 4:-1:1
        _get_col!(m, col, i)
        if sum(norm.(col)) > 0.5
            col_id = i
            break
        end
    end

    for i in 4:-1:1
        _get_row!(m, row, i)
        if sum(norm.(row)) > 0.5    
            row_id = i
            break
        end
    end

    col_val = sqrt(sum(col.^2))
    row_val = sqrt(sum(row.^2))

    for i in 0:3
        
        col_val = col_val * (-1)^(i)
        row_val = row_val * (-1)^(i÷2)

        a, b, c, d = col./col_val
        p, q, r, s = row./row_val

        if col[row_id]/col_val ≈ row_val && row[col_id]/row_val ≈ col_val
            @assert p^2 + q^2 + r^2 + s^2 ≈ 1
            @assert a^2 + b^2 + c^2 + d^2 ≈ 1
            A = [p + q * im s + r * im; -s + r * im p - q * im]
            B = [a + b * im -d + c * im; d + c * im a - b * im]
            return A, B    
        end
    end
    error("Algorithm does not work")
end

function decomposeSO4(m::AbstractMatrix{ComplexF64})
    decomposeSO4(real.(m))
end

"""
    compile2qubit(m::AbstractMatrix{ComplexF64}, N, wires=nothing)

compiles an arbitrary U(4) matrix into a quantum circuit. Algorithm taken from
arxiv:quant-ph/0308006, arxiv:quant-ph/0211002
"""
function compile2qubit(m::AbstractMatrix{<:Complex}, N, wires=nothing)
    cg = CircuitGate[]
    E = 1 / sqrt(2) .* [1 im 0 0; 0 0 im 1; 0 0 im -1; 1 -im 0 0]
    U = E' * m * E
    P2 = U * transpose(U)
    Diag, K_2 = eigen(P2)
    
    K_2 = real.(K_2)
    K_2 .= gramm_schmidt!(K_2)
    Diag[1] = Diag[1] / det(K_2)^2
    K_2[:,1] = K_2[:,1] * det(K_2)
    
    Diag = diagm(sqrt.(Diag))
    Diag[1] = Diag[1] * det(U) / det(Diag)

    P = K_2 * Diag * inv(K_2)
    K_1 = inv(P) * U
    C, D = decomposeSO4(inv(K_2) * K_1)

    # println(norm(E'*kron(C,D)*E - inv(K_2) * K_1))
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
                    circuit_gate(2, X, 1)])

    A, B = decomposeSO4(K_2)

    # println(norm(E'*kron(A,B)*E - K_2))
    append!(cg, compile1qubit(A, [2])[1])
    append!(cg, compile1qubit(B, [1])[1])

    return cg
end

"""
    fillψ!(d::Vector{ComplexF64}, ψ::Vector{Float64}, l::Int)

fills empty vector `ψ`, such that  ψ = −im*[logχ1(d) log χ2(d) ··· logχl−1(U)],
where χn(d) = d[2n-1]*d[2n+2]/(d[2n]*d[2n+1]) per arxiv:quant-ph/0303039
"""
function fillψ!(d::AbstractVector{<:Complex}, ψ::Vector{<:Real}, l::Int)
    for i in StepRange(1, 1, l - 1)
        ψ[i] = imag(log(d[2i - 1] * d[2i + 2] / (d[2i] * d[2i + 1])))
    end
end

"""
    ηcol(l::Int, n::Int)

calculates a single flip state for given Int `n` and a column length of `l`,
    where 0 <= n <= l. returns a vector of length l, where vec[n] = 1, vec[n+1]
    = -1. if n == 0, vec[1] = -1, if n==l, vec[l] = 1
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
@memoize Dict function ηcol(l::Int, n::Int)
    vec = zeros(Int, l)
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
    flip_state(m::Int, l::Int)

calculates the `m`th column of the η⊗ matrix by summing all flip states
associated with the binary representation of `m` over the span of [1,l-1]
algorithm taken from arxiv:quant-ph/0303039
"""
@memoize Dict function flip_state(m::Int, l::Int)
    # calculates flip states of a given Int `x`
    f = filter(x -> count_ones(x & m) % 2 == 1, 1:l - 1)
    sum(ηcol.((l - 1,), f))
end

"""
    svalue(i::Int)

returns tuple of position of 1s in the binary representation of Int `i`
starting from least significant digit.
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
function svalue(i::Int)
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

recursively compiles an arbitrary diagonal unitary matrix D_{2^N} ∈ C^(2^N×2^N)
into a quantum circuit. algorithm taken from arxiv:quant-ph/0303039 D_{2^N} can
be decomposed into D_{2^(N-1)} ∈ C^(2^(N-1)×2^(N-1)) ⊗ exp(i*ϕ)*RzGate(θ)
"""
function compilediag(d::Vector{<:Complex}, N, cg=nothing, j=0)
    cgs = CircuitGate[]
    # single qubit diagonal matrix can be split into a phase shift + RzGate
    if length(d) == 2
        logd = imag.(log.(d))
        β = (logd[1] - logd[2])
        if isnothing(cg)
            cg = CircuitGate[]
        end
        if !(β ≈ 0)
            push!(cg, circuit_gate((j + 1,), RzGate(-β)))
        end
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
        if !(α[i]≈0)
            push!(cgs, circuit_gate((j + 1,), RzGate(α[i])))
        end
        rz[s] .= 1 ./ rz[s]
        d = d ./ rz
    end

    for k in reverse(oldwires)
        push!(cgs, circuit_gate(j + 1, X, k + 1 + j))
    end

    if !(d[1] ≈ d[2])
        logd = imag.(log.(d))
        α[1] = logd[1] - logd[2]
        if !(α[1] ≈ 0)
            pushfirst!(cgs, circuit_gate((j + 1,), RzGate(-α[1])))
        end
        d[1:2:end] = d[1:2:end] ./ exp(im * α[1] / 2)
    end
    if isnothing(cg)
        return compilediag(d[1:2:end], N, cgs, j + 1)
    end
    return compilediag(d[1:2:end], N, append!(cg, cgs), j + 1)
end
