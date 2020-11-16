using Memoize

"""
    compile!(m::Matrix{Complex64,2}) where {N}

compiles quantum circuit from given unitary matrix.
"""

function compile(m::Matrix{Complex64,2}, N, wires=nothing)
    isunitary(m) || error("Only unitary matrices can be compiled into a quantum circuit")
    s = size(m)
    s[1] == s[2] || error("Only square matrices can be compiled into a quantum circuit")

    if Diagonal(diag(m)) ≈ m
        return compilediag(diag(m), N, wires)[1]
    end

    if size(m)[1] == 2
        return compile1qubit(m, N, wires)[1]
    end

    if size(m)[1] == 4
        return compile2qubit(m, N, wires)[1]
    end

end

function compile1qubit(m::Matrix{ComplexF64}, N, wires=nothing)
    phase = sqrt(m[1,1]*m[2,2]/norm(m[1,1]*m[2,2]))
    b = m / phase
    θ1 = imag(log(2im*b[2,1]*b[2,2]/sin(angle2)))
    ϕ = real(acos(sqrt(b[1,1]*b[2,2]))*2)
    θ2 = imag(-log(2*im*b[1,1]*b[2,1]/sin(angle2)))
    if isnothing(wires)
        wires = [1]
    end
    return CircuitGateChain{N}(single_qubit_circuit_gate.((wires[1],), [RzGate(θ2), RxGate(ϕ), RzGate(θ1)], (N,))), phase
end


function η!(d::Vector{ComplexF64}, ψ::Vector{Float64}, l::Integer)
    for i in StepRange(1,1,l-1)
        ψ[i] = imag(log(d[2i-1]*d[2i+2]/(d[2i]*d[2i+1])))
    end
end

grayencode(n::Integer) = n ⊻ (n >> 1)

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
    vec[n+1] = -1
    vec
end

@memoize Dict function flip_code(m::Integer, l::Integer)
    f = filter(x->count_ones(x&m)%2==1, 1:l-1)
    sum(ηcol.((l-1,), f))
end

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

#Algorithm taken from https://arxiv.org/pdf/quant-ph/0303039.pdf
function compilediag(d::Vector{ComplexF64}, N, cg=nothing, j=0)
    cgs = AbstractCircuitGate{N}[]
    if length(d) == 2
        logd = imag.(log.(d))
        β = (logd[1] - logd[2])
        if isnothing(cg)
            cg = AbstractCircuitGate{N}[]
        end
        push!(cg, CircuitGate((j+1,), RzGate(-β), N))
        return cg
    end
    l = 2^(N-j-1)
    ψ = zeros(Float64, l-1)
    α = zeros(Float64, l)

    η!(d, ψ, l)

    ηplus = zeros(Float64, (l-1, l-1))
    for i in StepRange(1,1,l-1)
        ηplus[:,i] = flip_code(grayencode(i), l)
    end

    ηplus = inv(ηplus)
    α[2:end] = -0.5 * ηplus * ψ

    gatewires = svalue(grayencode(1))
    oldwires = ()

    for i in StepRange(2, 1, l)
        rz = [exp(-im*α[i]/2), exp(im*α[i]/2)]
        for _ in 1:N-j-1
            rz = append!(rz, rz)
        end
        s = filter(x->count_ones(x&grayencode(i-1))%2==1, 1:l-1)
        s = append!(2 .* s .+ 1, 2 .* s .+ 2)

        for k in symdiff(gatewires, oldwires)
            push!(cgs, CircuitGate((j+1, k+1+j), ControlledGate{1,2}(X), N))
        end
        oldwires = gatewires
        gatewires = svalue(grayencode(i))
        push!(cgs, CircuitGate((j+1,), RzGate(α[i]), N))

        rz[s] .= 1 ./rz[s]
        d = d ./ rz
    end

    for k in reverse(oldwires)
        push!(cgs, CircuitGate((j+1, k+1+j), ControlledGate{1,2}(X), N))
    end

    if !(d[1] ≈ d[2])
        logd = imag.(log.(d))
        α[1] = logd[1] - logd[2]
        pushfirst!(cgs, CircuitGate((j+1,), RzGate(-α[1]), N))
        d[1:2:end] = d[1:2:end] ./ exp(im*α[1]/2)
    end
    if isnothing(cg)
        return compilediag(d[1:2:end], N, cgs, j+1)
    end
    return compilediag(d[1:2:end], N, append!(cg, cgs), j+1)
end
