using LinearAlgebra
using StaticArrays
using SparseArrays

"""
    AbstractGate

Abtract unitary quantum gate. `N` is the number of "wires" the gate acts on.
"""
abstract type AbstractGate end
matrix(g::AbstractGate)::SparseMatrixCSC{Complex{Float64},Int} = matrix(typeof(g), data(g))

"""
Pauli X gate

``X = \\begin{pmatrix} 0 & 1 \\\\ 1 & 0 \\end{pmatrix}``
"""
struct XGate <: AbstractGate end

"""
Pauli Y gate

``Y = \\begin{pmatrix} 0 & -i \\\\ i & 0 \\end{pmatrix}``
"""
struct YGate <: AbstractGate end

"""
Pauli Z gate

``Z = \\begin{pmatrix} 1 & 0 \\\\ 0 & -1 \\end{pmatrix}``
"""
struct ZGate <: AbstractGate end

matrix(::Type{XGate}, ::Nothing) = sparse([1,2],[2,1], ComplexF64[1,1])
matrix(::Type{YGate}, ::Nothing) = sparse([1,2],[2,1], ComplexF64[-im,im])
matrix(::Type{ZGate}, ::Nothing) = sparse([1,2],[1,2], ComplexF64[1,-1])

LinearAlgebra.ishermitian(::XGate) = true
LinearAlgebra.ishermitian(::YGate) = true
LinearAlgebra.ishermitian(::ZGate) = true

# Pauli matrices are Hermitian
Base.adjoint(::XGate) = X
Base.adjoint(::YGate) = Y
Base.adjoint(::ZGate) = Z

# corresponding instances
X = XGate()
Y = YGate()
Z = ZGate()

# wires
num_wires(::XGate)::Int = 1
num_wires(::YGate)::Int = 1
num_wires(::ZGate)::Int = 1


"""
Hadamard gate

``H = \\frac{\\sqrt} \\begin{pmatrix} 1 & 1 \\\\ 1 & 1 \\end{pmatrix}``
"""
struct HadamardGate <: AbstractGate end

matrix(::Type{HadamardGate}, ::Nothing) = ComplexF64[1 1; 1 -1] / sqrt(2)

LinearAlgebra.ishermitian(::HadamardGate) = true
# Hadamard gate is Hermitian
Base.adjoint(H::HadamardGate) = H

# wires
num_wires(::HadamardGate)::Int = 1


"""
S gate

``S = \\frac{\\sqrt} \\begin{pmatrix} 1 & 0 \\\\ 0 & i \\end{pmatrix}``
"""
struct SGate <: AbstractGate end


"""
T gate

``T = \\frac{\\sqrt} \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{\\frac{iπ}{4}} \\end{pmatrix}``
"""
struct TGate <: AbstractGate end


"""
S† gate

``S^{†} = \\begin{pmatrix} 1 & 0 \\\\ 0 & -i \\end{pmatrix}``
"""
struct SdagGate <: AbstractGate end

"""
T† gate

``T^{†} = \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{-\\frac{iπ}{4}} \\end{pmatrix}``
"""
struct TdagGate <: AbstractGate end

matrix(::Type{SGate}, ::Nothing)::Matrix{ComplexF64} = sparse([1,2],[1,2], ComplexF64[1,im])
matrix(::Type{TGate}, ::Nothing)::Matrix{ComplexF64} = sparse([1,2],[1,2], ComplexF64[1,Base.exp(im * π / 4)])

matrix(::Type{SdagGate}, ::Nothing)::Matrix{ComplexF64} = sparse([1,2],[1,2], ComplexF64[1,-im])
matrix(::Type{TdagGate}, ::Nothing)::Matrix{ComplexF64} = sparse([1,2],[1,2], ComplexF64[1,Base.exp(-im * π / 4)])

LinearAlgebra.ishermitian(::SGate) = false
LinearAlgebra.ishermitian(::TGate) = false

LinearAlgebra.ishermitian(::SdagGate) = false
LinearAlgebra.ishermitian(::TdagGate) = false

Base.adjoint(::SGate) = SdagGate()
Base.adjoint(::TGate) = TdagGate()

Base.adjoint(::SdagGate) = SGate()
Base.adjoint(::TdagGate) = TGate()

# wires
num_wires(::SGate)::Int = 1
num_wires(::TGate)::Int = 1
num_wires(::SdagGate)::Int = 1
num_wires(::TdagGate)::Int = 1


"""
Rotation-X gate

``R_{x}(\\theta) = \\begin{pmatrix} \\cos(\\frac{\\theta}) & -i\\sin(\\frac{\\theta}) \\\\ -i\\sin(\\frac{\\theta}) & \\cos(\\frac{\\theta}) \\end{pmatrix}``
"""
struct RxGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{Float64}

    function RxGate(θ::Real)
        new([θ])
    end
end


function matrix(::Type{RxGate}, data::Vector{Float64})
    c = cos(data[] / 2.0)
    s = sin(data[] / 2.0)
    sparse([1,2,1,2],[1,1,2,2], ComplexF64[c, -im * s, -im * s, c])
end



LinearAlgebra.ishermitian(g::RxGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::RxGate)::Int = 1


"""
Rotation-Y gate

``R_{y}(\\theta) = \\begin{pmatrix} \\cos(\\frac{\\theta}) & -\\sin(\\frac{\\theta}) \\\\ \\sin(\\frac{\\theta}) & \\cos(\\frac{\\theta}) \\end{pmatrix}``
"""
struct RyGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{Float64}

    function RyGate(θ::Real)
        new([θ])
    end
end

function matrix(::Type{RyGate}, data::Vector{Float64})
    c = cos(data[] / 2.0)
    s = sin(data[] / 2.0)
    sparse([1,2,1,2],[1,1,2,2], ComplexF64[c, s, -s, c])
end

LinearAlgebra.ishermitian(g::RyGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::RyGate)::Int = 1

"""
Rotation-Z gate

``R_{z}(\\theta) = \\begin{pmatrix} e^{\\frac{-i\\theta}} & 0 \\\\ 0 & e^{\\frac{i\\theta}} \\end{pmatrix}``
"""
struct RzGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{Float64}
    function RzGate(θ::Real)
        new([θ])
    end
end

function matrix(::Type{RzGate}, data::Vector{Float64})
    eθ = exp(im * data[] / 2)
    sparse([1,2],[1,2], ComplexF64[conj(eθ), eθ])
end

LinearAlgebra.ishermitian(g::RzGate) = abs(sin(g.θ[] / 2)) < 4 * eps()


Base.adjoint(g::RxGate) = RxGate(-g.θ[])
Base.adjoint(g::RyGate) = RyGate(-g.θ[])
Base.adjoint(g::RzGate) = RzGate(-g.θ[])

# wires
num_wires(::RzGate)::Int = 1

"""
General rotation operator gate: rotation by angle `θ` around unit vector `n`.

``R_{\\vec}(\\theta) = \\cos(\\frac{\\theta})I - i\\sin(\\frac{\\theta})\\vec\\sigma, \\\\ \\sigma = [X, Y, Z]``
"""
struct RotationGate <: AbstractGate
    nθ::Vector{Float64}

    function RotationGate(nθ::Vector{<:Real})
        length(nθ) == 3 || error("Rotation axis vector must have length 3.")
        new(nθ)
    end

    function RotationGate(θ::Real, n::Vector{<:Real})
        length(n) == 3 || error("Rotation axis vector must have length 3.")
        norm(n) ≈ 1 || error("Norm of rotation axis vector must be 1.")
        new(n * θ)
    end
end

function matrix(::Type{RotationGate}, data::Vector{Float64})
    θ = norm(data)
    if θ == 0
        return sparse([1,2],[1,2], ComplexF64[1, 1])
    end
    n = data / θ
    s = im * sin(θ / 2)
    c = cos(θ / 2)
    v = ComplexF64[c - s*n[3], - s*n[1] - im * s*n[2], - s*n[1] + im * s*n[2], cos(θ / 2) + s*n[3]]
    sparse([1,2,1,2],[1,1,2,2], v)
end

LinearAlgebra.ishermitian(g::RotationGate) = abs(sin(norm(g.nθ) / 2)) < 8 * eps()

Base.adjoint(g::RotationGate) = RotationGate(-g.nθ)

# wires
num_wires(::RotationGate)::Int = 1


"""
Phase shift gate

``P(\\phi) = \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{i\\phi} \\end{pmatrix}``
"""
struct PhaseShiftGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    ϕ::Vector{Float64}

    function PhaseShiftGate(ϕ::Real)
        new([ϕ])
    end
end

matrix(::Type{PhaseShiftGate}, data::Vector{Float64}) = sparse([1,2],[1,2], ComplexF64[1, Base.exp(im * data[])])

function LinearAlgebra.ishermitian(g::PhaseShiftGate)
    if abs(g.ϕ[]) < eps()
        return true
    end
    return false
end

Base.adjoint(g::PhaseShiftGate) = PhaseShiftGate(-g.ϕ[])

# wires
num_wires(::PhaseShiftGate)::Int = 1

"""
Swap gate

``SWAP = \\begin{pmatrix} 1 & 0 & 0 & 0 \\\\ 0 & 0 & 1 & 0 \\\\ 0 & 1 & 0 & 0 \\\\ 0 & 0 & 0 & 1 \\end{pmatrix}``
"""
struct SwapGate <: AbstractGate end


# matrix(::Type{SwapGate}, ::Nothing) = ComplexF64[1. 0. 0. 0.; 0. 0. 1. 0.; 0. 1. 0. 0.; 0. 0. 0. 1.]
matrix(::Type{SwapGate}, ::Nothing) = sparse([1,3,2,4],[1,2,3,4], ComplexF64[1,1,1,1])

# swap gate is Hermitian
LinearAlgebra.ishermitian(::SwapGate) = true
Base.adjoint(s::SwapGate) = s

# wires
num_wires(::SwapGate)::Int = 2

"""
Entanglement-XX gate

``G_{x}(\\theta) = e^{-i \\theta X \\otimes X / 2}``

Reference:\n
    B. Kraus and J. I. Cirac\n
    Optimal creation of entanglement using a two-qubit gate\n
    Phys. Rev. A 63, 062309 (2001)
"""
struct EntanglementXXGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{Float64}
    function EntanglementXXGate(θ::Real)
        new([θ])
    end
end

function matrix(::Type{EntanglementXXGate}, data::Vector{Float64})
    c = cos(data[] / 2)
    s = sin(data[] / 2)
    sparse([1,4,2,3,2,3,1,4],[1,1,2,2,3,3,4,4], ComplexF64[c,-im*s,c,-im*s,-im*s,c,-im*s,c])
end

LinearAlgebra.ishermitian(g::EntanglementXXGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::EntanglementXXGate)::Int = 2

"""
Entanglement-YY gate

``G_{y}(\\theta) = e^{-i \\theta Y \\otimes Y / 2}``

Reference:\n
    B. Kraus and J. I. Cirac\n
    Optimal creation of entanglement using a two-qubit gate\n
    Phys. Rev. A 63, 062309 (2001)
"""
struct EntanglementYYGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{Float64}
    function EntanglementYYGate(θ::Real)
        new([θ])
    end
end

function matrix(::Type{EntanglementYYGate}, data::Vector{Float64})
    c = cos(data[] / 2)
    s = sin(data[] / 2)
    sparse([1,4,2,3,2,3,1,4],[1,1,2,2,3,3,4,4], ComplexF64[c,im*s,c,-im*s,-im*s,c,im*s,c])
end

LinearAlgebra.ishermitian(g::EntanglementYYGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::EntanglementYYGate)::Int = 2

"""
Entanglement-ZZ gate

``G_{z}(\\theta) = e^{-i \\theta Z \\otimes Z / 2}``

Reference:\n
    B. Kraus and J. I. Cirac\n
    Optimal creation of entanglement using a two-qubit gate\n
    Phys. Rev. A 63, 062309 (2001)
"""
struct EntanglementZZGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{Float64}
    function EntanglementZZGate(θ::Real)
        new([θ])
    end
end

function matrix(::Type{EntanglementZZGate}, data::Vector{Float64})
    eθ = exp(im * data[] / 2)
    sparse([1,2,3,4],[1,2,3,4], ComplexF64[conj(eθ), eθ, eθ, conj(eθ)])
end


LinearAlgebra.ishermitian(g::EntanglementZZGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::EntanglementZZGate)::Int = 2

Base.adjoint(g::EntanglementXXGate) = EntanglementXXGate(-g.θ[])
Base.adjoint(g::EntanglementYYGate) = EntanglementYYGate(-g.θ[])
Base.adjoint(g::EntanglementZZGate) = EntanglementZZGate(-g.θ[])


"""
General controlled gate: the `M` wires corresponding to the fastest running indices are the target and the remaining `N - M` wires the control
"""
struct ControlledGate{G} <: AbstractGate
    U::G
    M::Int

    function ControlledGate(U::AbstractGate, M::Int)
        new{typeof(U)}(U, M)
    end
end

# wires
num_wires(g::ControlledGate)::Int = g.M + num_wires(g.U)
target_wires(g::ControlledGate) = num_wires(g.U)
control_wires(g::ControlledGate) = g.M

function matrix(g::ControlledGate{G}) where {G <:AbstractGate}
    N = num_wires(g)
    Umat = matrix(g.U)
    CU = sparse(1:2^N, 1:2^N, ones(ComplexF64, 2^N))
    # Note: target qubit(s) corresponds to fastest varying index
    CU[end - size(Umat, 1) + 1:end, end - size(Umat, 2) + 1:end] = Umat
    return dropzeros!(CU)
end

LinearAlgebra.ishermitian(g::ControlledGate) =
    LinearAlgebra.ishermitian(g.U)

Base.adjoint(g::ControlledGate) =
    ControlledGate(Base.adjoint(g.U), g.M)

controlled_not() = ControlledGate(X, 1)

isunitary(m::AbstractMatrix) = (m * Base.adjoint(m) ≈ I)


"""
MatrixGate: general gate constructed from an unitary matrix
"""
struct MatrixGate <: AbstractGate
    matrix::SparseMatrixCSC{Complex{Float64},Int}
    function MatrixGate(m)
        @assert size(m, 1) == size(m, 2)
        isunitary(m) || error("Quantum operators must be unitary")
        return new(sparse(m))
    end
end

# wires
num_wires(g::MatrixGate)::Int = Int(log(2, size(g.matrix, 1)))

matrix(g::MatrixGate) = g.matrix

matrix(g::SparseMatrixCSC{Complex{Float64},Int}) = g

Base.adjoint(g::MatrixGate) = MatrixGate(Base.adjoint(g.matrix))

LinearAlgebra.ishermitian(g::MatrixGate) =
    LinearAlgebra.ishermitian(g.matrix)

function Base.isapprox(g1::G, g2::G) where {G <: AbstractGate} 
    for name in fieldnames(G)
        if !(getfield(g1, name) ≈ getfield(g2, name))
            return false
        end
    end
    return true
end

# handle different gate types or dimensions
Base.isapprox(::AbstractGate, ::AbstractGate) = false

data(g::AbstractGate)::Union{Vector{Float64},Nothing} = get_data(g)
get_data(::XGate) = nothing
get_data(::YGate) = nothing
get_data(::ZGate) = nothing
get_data(::HadamardGate) = nothing
get_data(::SGate) = nothing
get_data(::SdagGate) = nothing
get_data(::TGate) = nothing
get_data(::TdagGate) = nothing
get_data(::SwapGate) = nothing
get_data(g::RxGate) = g.θ
get_data(g::RyGate) = g.θ
get_data(g::RzGate) = g.θ
get_data(g::RotationGate) = g.nθ
get_data(::MatrixGate) = nothing
get_data(g::PhaseShiftGate) = g.ϕ
get_data(g::ControlledGate) = data(g.U)
get_data(g::EntanglementXXGate) = g.θ
get_data(g::EntanglementYYGate) = g.θ
get_data(g::EntanglementZZGate) = g.θ
