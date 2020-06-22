using LinearAlgebra

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

"""
Pauli X Matrix

``X = \\begin{pmatrix} 0 & 1 \\\\ 1 & 0 \\end{pmatrix}``
"""
struct XGate <: AbstractGate{1} end
"""
Pauli Y Matrix

``Y = \\begin{pmatrix} 0 & -i \\\\ i & 0 \\end{pmatrix}``
"""
struct YGate <: AbstractGate{1} end

"""
Pauli Z Matrix

``Z = \\begin{pmatrix} 1 & 0 \\\\ 0 & -1 \\end{pmatrix}``
"""
struct ZGate <: AbstractGate{1} end

matrix(::XGate) = [0.  1.; 1.  0.]
matrix(::YGate) = [0. -im; im  0.]
matrix(::ZGate) = [1.  0.; 0. -1.]

LinearAlgebra.ishermitian(::XGate) = true
LinearAlgebra.ishermitian(::YGate) = true
LinearAlgebra.ishermitian(::ZGate) = true

# Pauli matrices are Hermitian
Base.adjoint(X::XGate) = X
Base.adjoint(Y::YGate) = Y
Base.adjoint(Z::ZGate) = Z

# corresponding instances
X = XGate()
Y = YGate()
Z = ZGate()


"""
Hadamard Matrix

``H = \\frac{1}{\\sqrt{2}} \\begin{pmatrix} 1 & 1 \\\\ 1 & 1 \\end{pmatrix}``
"""
struct HadamardGate <: AbstractGate{1} end

matrix(::HadamardGate) = [1 1; 1 -1] / sqrt(2)

LinearAlgebra.ishermitian(::HadamardGate) = true
# Hadamard gate is Hermitian
Base.adjoint(H::HadamardGate) = H


# S & T gates
"""
S Matrix

``S = \\frac{1}{\\sqrt{2}} \\begin{pmatrix} 1 & 0 \\\\ 0 & i \\end{pmatrix}``
"""
struct SGate <: AbstractGate{1} end
"""
T Matrix

``T = \\frac{1}{\\sqrt{2}} \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{\\frac{iπ}{4}} \\end{pmatrix}``
"""
struct TGate <: AbstractGate{1} end

"""
S† Matrix

``S^{†} = \\begin{pmatrix} 1 & 0 \\\\ 0 & -i \\end{pmatrix}``
"""
struct SdagGate <: AbstractGate{1} end
"""
T† Matrix

``T^{†} = \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{-\\frac{iπ}{4}} \\end{pmatrix}``
"""
struct TdagGate <: AbstractGate{1} end

matrix(::SGate) = [1. 0.; 0. im]
matrix(::TGate) = [1. 0.; 0. exp(im*π/4)]

matrix(::SdagGate) = [1. 0.; 0. -im]
matrix(::TdagGate) = [1. 0.; 0. exp(-im*π/4)]

LinearAlgebra.ishermitian(::SGate) = false
LinearAlgebra.ishermitian(::TGate) = false

LinearAlgebra.ishermitian(::SdagGate) = false
LinearAlgebra.ishermitian(::TdagGate) = false

Base.adjoint(::SGate) = SdagGate()
Base.adjoint(::TGate) = TdagGate()

Base.adjoint(::SdagGate) = SGate()
Base.adjoint(::TdagGate) = TGate()


# rotation gates
"""
Rotation X Matrix

``R_{x}(\\theta) = \\begin{pmatrix} \\cos(\\frac{\\theta}{2}) & -i\\sin(\\frac{\\theta}{2}) \\\\ -i\\sin(\\frac{\\theta}{2}) & \\cos(\\frac{\\theta}{2}) \\end{pmatrix}``
"""
struct RxGate <: AbstractGate{1}
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{<:Real}

    function RxGate(θ::Real)
        new([θ])
    end
end

function matrix(g::RxGate)
    c = cos(g.θ[]/2)
    s = sin(g.θ[]/2)
    [c -im*s; -im*s c]
end

function LinearAlgebra.ishermitian(g::RxGate)
    if mod2pi(g.θ[]) < eps()
        return true
    end
    return false
end

"""
Rotation Y Matrix

``R_{y}(\\theta) = \\begin{pmatrix} \\cos(\\frac{\\theta}{2}) & -\\sin(\\frac{\\theta}{2}) \\\\ \\sin(\\frac{\\theta}{2}) & \\cos(\\frac{\\theta}{2}) \\end{pmatrix}``
"""
struct RyGate <: AbstractGate{1}
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{<:Real}

    function RyGate(θ::Real)
        new([θ])
    end
end

function matrix(g::RyGate)
    c = cos(g.θ[]/2)
    s = sin(g.θ[]/2)
    [c -s; s c]
end

function LinearAlgebra.ishermitian(g::RyGate)
    if mod2pi(g.θ[]) < eps()
        return true
    end
    return false
end

"""
Rotation Z Matrix

``R_{z}(\\theta) = \\begin{pmatrix} e^{\\frac{-i\\theta}{2}} & 0 \\\\ 0 & e^{\\frac{i\\theta}{2}} \\end{pmatrix}``
"""
struct RzGate <: AbstractGate{1}
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{<:Real}
    function RzGate(θ::Real)
        new([θ])
    end
end

function matrix(g::RzGate)
    [exp(-im*g.θ[]/2) 0; 0 exp(im*g.θ[]/2)]
end

function LinearAlgebra.ishermitian(g::RzGate)
    if mod2pi(g.θ[]) < eps()
        return true
    end
    return false
end

Base.adjoint(g::RxGate) = RxGate(-g.θ[])
Base.adjoint(g::RyGate) = RyGate(-g.θ[])
Base.adjoint(g::RzGate) = RzGate(-g.θ[])

# general rotation operator gate
"""
General Rotation Matrix
Rotation by angle `θ` around unit vector `n⃗`.

``R_{\\vec{n}}(\\theta) = \\cos(\\frac{\\theta}{2})I - i\\sin(\\frac{\\theta}{2})\\vec{n}\\sigma, \\\\ \\sigma = [X, Y, Z]``
"""
struct RotationGate <: AbstractGate{1}
    nθ::AbstractVector{<:Real}

    function RotationGate(nθ::AbstractVector{<:Real})
        length(nθ) == 3 || error("Rotation axis vector must have length 3.")
        new(nθ)
    end

    function RotationGate(θ::Real, n::AbstractVector{<:Real})
        length(n) == 3 || error("Rotation axis vector must have length 3.")
        norm(n) ≈ 1 || error("Norm of rotation axis vector must be 1.")
        new(n*θ)
    end
end

function matrix(g::RotationGate)
    θ = norm(g.nθ)
    if θ == 0
        return Matrix{Complex{eltype(g.nθ)}}(I, 2, 2)
    end
    n = g.nθ/θ
    cos(θ/2)*I - im*sin(θ/2)*pauli_vector(n...)
end

function LinearAlgebra.ishermitian(g::RotationGate)
    if norm(g.nθ + g.nθ) < eps()
        return true
    end
    return false
end

Base.adjoint(g::RotationGate) = RotationGate(-g.nθ)


"""
Phase Shift Gate

``P(\\phi) = \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{i\\phi} \\end{pmatrix}``
"""
struct PhaseShiftGate <: AbstractGate{1}
    # use a reference type (array with 1 entry) for compatibility with Flux
    ϕ::Vector{<:Real}

    function PhaseShiftGate(ϕ::Real)
        new([ϕ])
    end
end

matrix(g::PhaseShiftGate) = [1 0; 0 exp(im*g.ϕ[])]

function LinearAlgebra.ishermitian(g::PhaseShiftGate)
    if g.ϕ[] == 0
        return true
    end
    return false
end

Base.adjoint(g::PhaseShiftGate) = PhaseShiftGate(-g.ϕ[])



# swap gate
"""
Swap Gate

``SWAP = \\begin{pmatrix} 1 & 0 & 0 & 0 \\\\ 0 & 0 & 1 & 0 \\\\ 0 & 1 & 0 & 0 \\\\ 0 & 0 & 0 & 1 \\end{pmatrix}``
"""
struct SwapGate <: AbstractGate{2} end


matrix(::SwapGate) = [1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.]

# swap gate is Hermitian
LinearAlgebra.ishermitian(::SwapGate) = true
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

function LinearAlgebra.ishermitian(g::ControlledGate{M,N}) where {M,N}
    return LinearAlgebra.ishermitian(g.U)
end
Base.adjoint(g::ControlledGate{M,N}) where {M,N} = ControlledGate{M,N}(Base.adjoint(g.U))

controlled_not() = ControlledGate{1,2}(X)

struct MatrixGate{N} <: AbstractGate{N}
    matrix::AbstractMatrix
    function MatrixGate(m)
        d = 2
        @assert size(m,1) == size(m,2)
        N = Int(log(d, size(m,1)))
        return new{N}(m)
    end
end

function matrix(MG::MatrixGate{N}) where N
    MG.matrix
end

function Base.adjoint(MG::MatrixGate{N}) where N
    return MatrixGate(Base.adjoint(MG.matrix))
end
