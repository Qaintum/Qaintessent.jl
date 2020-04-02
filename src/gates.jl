
"""
    pauli_vector(x, y, z)

Assemble the "Pauli vector" matrix.
"""
pauli_vector(x, y, z) = [z x-im*y; x+im*y -z]


"""
    AbstractGate{N}

Abtract unitary quantum gate. `N` is the number of "wires" the gate acts on.
"""
abstract type AbstractGate{N} end


# Pauli matrices

struct XGate <: AbstractGate{1} end
struct YGate <: AbstractGate{1} end
struct ZGate <: AbstractGate{1} end

matrix(::XGate) = [0.  1.; 1.  0.]
matrix(::YGate) = [0. -im; im  0.]
matrix(::ZGate) = [1.  0.; 0. -1.]

# Pauli matrices are Hermitian
Base.adjoint(X::XGate) = X
Base.adjoint(Y::YGate) = Y
Base.adjoint(Z::ZGate) = Z

# corresponding instances
X = XGate()
Y = YGate()
Z = ZGate()


# Hadamard gate

struct HadamardGate <: AbstractGate{1} end

matrix(::HadamardGate) = [1 1; 1 -1] / sqrt(2)

# Hadamard gate is Hermitian
Base.adjoint(H::HadamardGate) = H


# S & T gates

struct SGate <: AbstractGate{1} end
struct TGate <: AbstractGate{1} end

struct SdagGate <: AbstractGate{1} end
struct TdagGate <: AbstractGate{1} end

matrix(::SGate) = [1. 0.; 0. im]
matrix(::TGate) = [1. 0.; 0. exp(im*π/4)]

matrix(::SdagGate) = [1. 0.; 0. -im]
matrix(::TdagGate) = [1. 0.; 0. exp(-im*π/4)]

Base.adjoint(::SGate) = SdagGate()
Base.adjoint(::TGate) = TdagGate()

Base.adjoint(::SdagGate) = SGate()
Base.adjoint(::TdagGate) = TGate()


# rotation gates

struct RxGate <: AbstractGate{1}
    θ::Real
end

function matrix(g::RxGate)
    c = cos(g.θ/2)
    s = sin(g.θ/2)
    [c -im*s; -im*s c]
end

struct RyGate <: AbstractGate{1}
    θ::Real
end

function matrix(g::RyGate)
    c = cos(g.θ/2)
    s = sin(g.θ/2)
    [c -s; s c]
end

struct RzGate <: AbstractGate{1}
    θ::Real
end

function matrix(g::RzGate)
    [exp(-im*g.θ/2) 0; 0 exp(im*g.θ/2)]
end

Base.adjoint(g::RxGate) = RxGate(-g.θ)
Base.adjoint(g::RyGate) = RyGate(-g.θ)
Base.adjoint(g::RzGate) = RzGate(-g.θ)

# general rotation operator gate
struct RotationGate <: AbstractGate{1}
    nθ::AbstractVector{<:Real}

    function RotationGate(θ::Real, n::AbstractVector{<:Real})
        length(n) == 3 || error("Rotation axis vector must have length 3.")
        norm(n) ≈ 1 || error("Norm of rotation axis vector must be 1.")
        new(θ*n)
    end
end

function matrix(g::RotationGate)
    θ = norm(g.nθ)
    n = g.nθ/θ
    cos(θ/2)*I - im*sin(θ/2)*pauli_vector(n...)
end

Base.adjoint(g::RotationGate) = RotationGate(-norm(g.nθ), g.nθ/norm(g.nθ))


# phase shift gate

struct PhaseShiftGate <: AbstractGate{1}
    ϕ::Real
end

matrix(g::PhaseShiftGate) = [1 0; 0 exp(im*g.ϕ)]

Base.adjoint(g::PhaseShiftGate) = PhaseShiftGate(-g.ϕ)


# swap gate

struct SwapGate <: AbstractGate{2} end

matrix(::SwapGate) = [1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.]

# swap gate is Hermitian
Base.adjoint(s::SwapGate) = s


# general controlled gate

struct ControlledGate{M,N} <: AbstractGate{N}
    U::AbstractGate{M}
    function ControlledGate{M,N}(U::AbstractGate{M}) where {M,N}
        M < N || error("Number of target wires of a controlled gate must be smaller than overall number of wires.")
        new{M,N}(U)
    end
end

function matrix(g::ControlledGate{M,N}) where {M,N}
    Umat = matrix(g.U)
    CU = sparse(one(eltype(Umat))*I, 2^N, 2^N)
    # Note: following the ordering convention of `kron` here, i.e.,
    # second (target) qubit corresponds to fastest varying index
    CU[end-size(Umat,1)+1:end, end-size(Umat,2)+1:end] = Umat
    return CU
end

Base.adjoint(g::ControlledGate{M,N}) where {M,N} = ControlledGate{M,N}(Base.adjoint(g.U))

controlled_not() = ControlledGate{1,2}(X)
