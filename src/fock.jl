using LinearAlgebra
using QuantumOptics

"""
    FockState{M,N}

Quantum state represented as an `M` mode Fock basis truncated to `N` Fock states.
"""
struct FockState{M,N}
    "coefficients for various photon numbers in fock basis"
    state::AbstractArray{ComplexF64}
    "true if FockState is pure"
    pure::Bool

    function FockState{M,N}(state::AbstractArray{ComplexF64}, pure::Bool) where {M,N}
        size(state) == (M, N+1) || error("Dimension of `state` array must be (" * string(M) * ", " * string(N+1) * ")")
        new{M,N}(state, pure)
    end

    @doc """
        FockState(s::Int, m::Int, trunc::Int)

    Uniform quantum state represented as an `m` mode Fock basis truncated to `trunc` Fock states.
    """
    function FockState(s::Int, m::Int, trunc::Int)
        m > 0 || error("Number of modes in Fock state must be positive")
        trunc > 0 || error("Truncated photon number must be positive")
        trunc ≥ s ≥ 0 || error("Photon number must be positive and less than truncated photon number")

        state = zeros(ComplexF64, (m, trunc+1))
        state[:,s+1] .= 1.0 + 0.0im

        new{m,trunc}(state, true)
    end

    @doc """
        FockState(s::AbstractVector{Int}, m::Int, trunc::Int)

    General quantum state represented as an `m` mode Fock basis truncated to `trunc` Fock states.
    """
    function FockState(s::AbstractVector{Int}, m::Int, trunc::Int)
        m > 0 || error("Number of modes in Fock state must be positive")
        trunc > 0 || error("Truncated photon number must be positive")
        length(s) == m || error(string(m) * "-mode FockState created but " * string(length(s)) * " modes provided in `s`")
        all(trunc+1 .≥ s .≥ 0) || error("Photon number must be positive and less than truncated photon number")
        state = zeros(ComplexF64, (m, trunc+1))
        for i in 1:m
            state[i,s[i]+1] = 1.0 + 0.0im
        end

        new{m,trunc}(state, true)
    end
end

"""
    apply(cg::CircuitGate{K,N,G}, fs::FockState{M,N}; modes::Union{Nothing,AbstractVector{Int}}=nothing) where {K,M,N,G}

apply CircuitGate `cg` to FockState `fs`.
"""
function apply(cg::CircuitGate{N,N,G}, fs::FockState{M,N}; modes::Union{Nothing,AbstractVector{Int}}=nothing, normalize=true) where {M,N,G}
    if isnothing(modes)
        fs.state .= apply(Qaintessent.matrix(cg), fs.state; normalize=normalize)
        return fs
    end
    all(modes .> 0) || error("Selected modes `modes` must be greater than 0")
    fs.state[modes, :] .= apply(Qaintessent.matrix(cg), fs.state[modes, :]; normalize=normalize)
    fs
end

function apply(g::AbstractGate{N}, fs::FockState{M,N}; modes::Union{Nothing,AbstractVector{Int}}=nothing, normalize=true) where {M,N}
    if isnothing(modes)
        fs.state .= apply(Qaintessent.matrix(g), fs.state; normalize=normalize)
        return fs
    end
    all(modes .> 0) || error("Selected modes `modes` must be greater than 0")
    fs.state[modes, :] .= apply(Qaintessent.matrix(g), fs.state[modes, :]; normalize=normalize)
    fs
end

function apply(m::AbstractMatrix, fs_state::AbstractArray{<:Complex}; normalize=true)
    state = transpose(m*transpose(fs_state))
    if normalize
        for i in 1:size(state, 1)
            state[i,:] = state[i,:]/norm(state[i,:])
        end
    end
    state
end


function Base.isapprox(fs1::FockState{N, M}, fs2::FockState{N, M}) where {N,M}
    fs1.pure == fs2.pure && all(fs1.state == fs2.state)
end


function vacuum_fock_state(m::Int, trunc::Int)
    FockState(0, m, trunc)
end

function create_coherent_state(α::ComplexF64, trunc::Int)
    state = zeros(ComplexF64, trunc+1)
    state[1] = exp(-abs2(α)/2)
    for i in 1:trunc
        state[i+1] = state[i]*α/√i
    end
    state
end

function coherent_state(α::ComplexF64, m::Int, trunc::Int)
    state = transpose(create_coherent_state(α, trunc))
    FockState{m,trunc}(repeat(state; outer=[m,1]), true)
end

function coherent_state(αs::AbstractVector{ComplexF64}, m::Int, trunc::Int)
    length(αs) == m || error(string(m) * "-mode coherent state created but " * string(length(αs)) * " displacements `αs` provided")
    fs = vacuum_fock_state(m, trunc)
    state = zeros(ComplexF64, (m, trunc+1))
    for i in αs
        state[i,:] .= create_coherent_state(i, trunc)
    end
    FockState{m,trunc}(state, true)
end

coherent_state(α::Float64, trunc::Int) = coherent_state(complex(α, 0), 1, trunc)

function matrix(fs::FockState{N, M}) where {N,M}
    s = reduce(vcat, fs.state)
    kron(s, s)
end


"""
    â{N}

annihilation operator used in a Fock basis
"""

struct â{N} <:AbstractGate{N} end
Base.adjoint(::â{N})  where {N} = âDag{N}()

function matrix(annop::â{N}) where {N}
    m = zeros(ComplexF64,(N+1,N+1))
    m[diagind(m, 1)] .= sqrt.(1:N)
    return m
end

matrix(::CircuitGate{N,N,â{N}}) where {N} = matrix(â{N}())

"""
    âDag{N}

creation operator used in a Fock basis
"""

struct âDag{N} <:AbstractGate{N} end
Base.adjoint(::âDag{N}) where {N} = â{N}()

function matrix(creationop::âDag{N}) where {N}
    m = zeros(ComplexF64,(N+1,N+1))
    m[diagind(m, -1)] .= sqrt.(1:N)
    return m
end

matrix(::CircuitGate{N,N,âDag{N}}) where {N} = matrix(âDag{N}())


"""
    n̂{N}

number operator used in a Fock basis
"""
struct n̂{N} <:AbstractGate{N} end
Base.adjoint(::n̂{N}) where {N} = n̂{N}()

function matrix(numberop::n̂{N}) where {N}
    diagm(0:N)
end

matrix(::CircuitGate{N,N,n̂{N}}) where {N} = matrix(n̂{N}())

"""
    Û{N}

phase-shift operator used in a Fock basis
"""

struct Û{N} <:AbstractGate{N}
    θ::ComplexF64
end
Base.adjoint(u::Û{N}) where {N} = Û{N}(-u.θ)

function matrix(phaseop::Û{N}) where {N}
    exp(-phaseop.θ*im * Qaintessent.matrix(n̂{N}()))
end

matrix(cg::CircuitGate{N,N,Û{N}}) where {N} = matrix(Û{N}(cg.θ))

"""
    D̂{N}

displacement operator used in a Fock basis
"""

struct D̂{N} <:AbstractGate{N}
    α::ComplexF64
end
Base.adjoint(u::D̂{N}) where {N}  = D̂{N}(-u.α)

function matrix(displacementop::D̂{N}) where {N}
    exp(displacementop.α*Qaintessent.matrix(âDag{N}()) - conj(displacementop.α)*Qaintessent.matrix(â{N}()))
end

matrix(cg::CircuitGate{N,N,D̂{N}}) where {N} = matrix(D̂{N}(cg.α))

"""
    Ŝ{N}(ζ::Real)

squeezing operator used in a Fock basis
"""

struct Ŝ{N} <:AbstractGate{N}
    ζ::Real
end
Base.adjoint(s::Ŝ{N}) where {N}  = Ŝ{N}(-s.ζ)

function matrix(squeeze::Ŝ{N}) where {N}
    a2 = Qaintessent.matrix(â{N}())^2
    adag2 = Qaintessent.matrix(âDag{N}())^2
    exp(0.5squeeze.ζ * (a2 - adag2))
end

matrix(cg::CircuitGate{N,N,Ŝ{N}}) where {N} = matrix(Ŝ{N}(cg.ζ))

"""
    wigner()

calculates wigner function for given `xvec`, `yvec` using Clenshaw summation. algorithm taken from https://qojulia.org/
"""

function wigner(fs::FockState{M,N}, xvec::AbstractVector{<:Real}, yvec::AbstractVector{<:Real}; g=√2) where {M,N}
    statevector = vec(fs.state)
    rho = kron(statevector, adjoint(statevector))
    _2α = [complex(x, y)*sqrt(2) for x=xvec, y=yvec]
    abs2_2α = abs2.(_2α)
    w = zero(_2α)
    b0 = similar(_2α)
    b1 = similar(_2α)
    b2 = similar(_2α)

    @inbounds for L=N:-1:1
        QuantumOptics._clenshaw_grid(L, rho, abs2_2α, _2α, w, b0, b1, b2, 2)
    end

    QuantumOptics._clenshaw_grid(0, rho, abs2_2α, _2α, w, b0, b1, b2, 1)

    @inbounds for i=eachindex(w)
        abs2_2α[i] = exp(-abs2_2α[i]/2)/pi.*real(w[i])
    end
    abs2_2α
end
