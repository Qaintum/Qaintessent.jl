using LinearAlgebra

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
        new{M,N}(state, pure)
    end

    @doc """
        FockState(s::Int, m::Int, trunc::Int)

    Uniform quantum state represented as an `m` mode Fock basis truncated to `trunc` Fock states.
    """
    function FockState(s::Int, m::Int, trunc::Int)
        m > 0 || error("Number of modes in Fock state must be positive")
        trunc > 0 || error("Truncated photon number must be positive")
        trunc ≥ s > 0 || error("Photon number must be positive and less than truncated photon number")

        state = zeros(ComplexF64, (m, trunc))
        state[:,s] .= 1.0 + 0.0im

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
        all(trunc .≥ s .> 0) || error("Photon number must be positive and less than truncated photon number")
        state = zeros(ComplexF64, (m, trunc))
        for i in 1:m
            state[i,s[i]] = 1.0 + 0.0im
        end

        new{m,trunc}(state, true)
    end
end

"""
    apply(cg::CircuitGate{K,N,G}, fs::FockState{M,N}; modes::Union{Nothing,AbstractVector{Int}}=nothing) where {K,M,N,G}

apply CircuitGate `cg` to FockState `fs`.
"""
function apply(cg::CircuitGate{N,N,G}, fs::FockState{M,N}; modes::Union{Nothing,AbstractVector{Int}}=nothing) where {M,N,G}
    if isnothing(modes)
        fs.state .= apply(Qaintessent.matrix(cg), fs.state)
        return fs
    end
    all(modes .> 0) || error("Selected modes `modes` must be greater than 0")
    fs.state[modes, :] .= apply((Qaintessent.matrix(cg), fs.state[modes, :]))
    fs
end

function apply(g::AbstractGate{N}, fs::FockState{M,N}; modes::Union{Nothing,AbstractVector{Int}}=nothing) where {M,N}
    if isnothing(modes)
        fs.state .= apply(Qaintessent.matrix(g), fs.state)
        return fs
    end
    all(modes .> 0) || error("Selected modes `modes` must be greater than 0")
    fs.state[modes, :] .= apply(Qaintessent.matrix(g), fs.state[modes, :])
    fs
end

function apply(m::AbstractMatrix, fs_state::AbstractArray{<:Complex})
    transpose(m*transpose(fs_state))
end


function Base.isapprox(fs1::FockState{N, M}, fs2::FockState{N, M}) where {N,M}
    fs1.pure == fs2.pure && all(fs1.state == fs2.state)
end


function vacuum_fock_state(m::Int, trunc::Int)
    FockState(1, m, trunc)
end

function create_coherent_state(α::ComplexF64, trunc::Int)
    state = zeros(ComplexF64, trunc)
    state[1] = 1
    for i in 1:trunc-1
        state[i+1] = state[i]*α/√i
    end
    state /= norm(state)
end

function coherent_state(α::ComplexF64, m::Int, trunc::Int)
    state = create_coherent_state(α, trunc)
    FockState{m,trunc}([deepcopy(state) for i in 1:m], true)
end

function coherent_state(αs::AbstractVector{ComplexF64}, m::Int, trunc::Int)
    length(αs) == m || error(string(m) * "-mode coherent state created but " * string(length(αs)) * " displacements `αs` provided")
    fs = FockState(0, m, trunc)
    state = []
    for i in αs
        push!(state, create_coherent_state(i, trunc))
    end
    FockState{m,trunc}(state, true)
end

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
    m = zeros(ComplexF64,(N,N))
    m[diagind(m, 1)] .= sqrt.(1:N-1)
    return m
end


"""
    âDag{N}

creation operator used in a Fock basis
"""

struct âDag{N} <:AbstractGate{N} end
Base.adjoint(::âDag{N}) where {N} = â{N}()

function matrix(creationop::âDag{N}) where {N}
    m = zeros(ComplexF64,(N,N))
    m[diagind(m, -1)] .= sqrt.(1:N-1)
    return m
end

"""
    n̂{N}

number operator used in a Fock basis
"""
struct n̂{N} <:AbstractGate{N} end
Base.adjoint(::n̂{N}) where {N} = n̂{N}()

function matrix(numberop::n̂{N}) where {N}
    diagm(0:N-1)
end

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

"""
    D̂{N}

displacement operator used in a Fock basis
"""

struct D̂{N} <:AbstractGate{N}
    α::ComplexF64
end
Base.adjoint(u::D̂{N}) where {N}  = D̂{N}(-u.α)

function matrix(displacementop::D̂{N}) where {N}
    exp(displacementop.α*Qaintessent.matrix(âDag{N+1}()) - conj(displacementop.α)*Qaintessent.matrix(â{N+1}()))[1:N, 1:N]
end
