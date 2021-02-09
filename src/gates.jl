
"""
    AbstractGate

Abtract unitary quantum gate. `N` is the number of "wires" the gate acts on.
"""
abstract type AbstractGate end

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

matrix(::XGate) = ComplexF64[0  1 ; 1  0]
matrix(::YGate) = ComplexF64[0 -im; im 0]
matrix(::ZGate) = ComplexF64[1  0 ; 0 -1]

sparse_matrix(::XGate) = sparse(matrix(XGate()))
sparse_matrix(::YGate) = sparse(matrix(YGate()))
sparse_matrix(::ZGate) = sparse(matrix(ZGate()))

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

sparse_matrix(::HadamardGate) = sparse([1,1,2,2],[1,2,1,2], ComplexF64[1, 1, 1, -1] / sqrt(2)) 
matrix(::HadamardGate) = ComplexF64[1  1; 1 -1] / sqrt(2)

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

matrix(::SGate)    = ComplexF64[1 0; 0 im]
matrix(::SdagGate) = ComplexF64[1 0; 0 -im]
matrix(::TGate)    = ComplexF64[1 0; 0 Base.exp(im * π / 4)]
matrix(::TdagGate) = ComplexF64[1 0; 0 Base.exp(-im * π / 4)]

sparse_matrix(::SGate)    = sparse(matrix(SGate()))
sparse_matrix(::TGate)    = sparse(matrix(TGate()))
sparse_matrix(::SdagGate) = sparse(matrix(SdagGate()))
sparse_matrix(::TdagGate) = sparse(matrix(TdagGate()))

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

function matrix(g::RxGate)
    c = cos(g.θ[] / 2.0)
    s = sin(g.θ[] / 2.0)
    ComplexF64[c -im * s; -im * s c]
end

sparse_matrix(g::RxGate) = sparse(matrix(g))

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

function matrix(g::RyGate)
    c = cos(g.θ[] / 2.0)
    s = sin(g.θ[] / 2.0)
    ComplexF64[c -s; s c]
end

sparse_matrix(g::RyGate) = sparse(matrix(g))

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

function matrix(g::RzGate)
    eθ = exp(im * g.θ[] / 2)
    ComplexF64[conj(eθ) 0; 0 eθ]
end

sparse_matrix(g::RzGate) = sparse(matrix(g))

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

function matrix(g::RotationGate)
    θ = norm(g.nθ)
    if θ == 0
        return Matrix{ComplexF64}(I, 2, 2)
    end
    n = g.nθ / θ
    return cos(θ/2)*I - im*sin(θ/2)*pauli_vector(n...)
end

sparse_matrix(g::RotationGate) = sparse(matrix(g))

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

matrix(g::PhaseShiftGate) = ComplexF64[1 0; 0 Base.exp(im * g.ϕ[])]

sparse_matrix(g::PhaseShiftGate) = sparse(matrix(g))

LinearAlgebra.ishermitian(g::PhaseShiftGate) = abs(sin(g.ϕ[])) < 4 * eps()

Base.adjoint(g::PhaseShiftGate) = PhaseShiftGate(-g.ϕ[])

# wires
num_wires(::PhaseShiftGate)::Int = 1


"""
Swap gate

``SWAP = \\begin{pmatrix} 1 & 0 & 0 & 0 \\\\ 0 & 0 & 1 & 0 \\\\ 0 & 1 & 0 & 0 \\\\ 0 & 0 & 0 & 1 \\end{pmatrix}``
"""
struct SwapGate <: AbstractGate end

matrix(::SwapGate,) = ComplexF64[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

sparse_matrix(::SwapGate,) = sparse([1, 3, 2, 4], [1, 2, 3, 4], ComplexF64[1, 1, 1, 1])

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

function matrix(g::EntanglementXXGate)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    ComplexF64[c  0  0  -im*s; 0  c  -im*s  0; 0  -im*s  c  0; -im*s  0  0  c]
end

function sparse_matrix(g::EntanglementXXGate)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    sparse([1, 4, 2, 3, 2, 3, 1, 4], [1, 1, 2, 2, 3, 3, 4, 4], ComplexF64[c, -im*s, c, -im*s, -im*s, c, -im*s, c])
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

function matrix(g::EntanglementYYGate)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    ComplexF64[c  0  0  im*s; 0  c  -im*s  0; 0  -im*s  c  0; im*s  0  0  c]
end

function sparse_matrix(g::EntanglementYYGate)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    sparse([1, 4, 2, 3, 2, 3, 1, 4], [1, 1, 2, 2, 3, 3, 4, 4], ComplexF64[c, im*s, c, -im*s, -im*s, c, im*s, c])
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

function matrix(g::EntanglementZZGate)
    eθ = exp(im * g.θ[] / 2)
    ComplexF64[conj(eθ) 0 0 0; 0 eθ 0 0; 0 0 eθ 0; 0 0 0 conj(eθ)]
end

function sparse_matrix(g::EntanglementZZGate)
    eθ = exp(im * g.θ[] / 2)
    sparse([1, 2, 3, 4], [1, 2, 3, 4], ComplexF64[conj(eθ), eθ, eθ, conj(eθ)])
end

LinearAlgebra.ishermitian(g::EntanglementZZGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::EntanglementZZGate)::Int = 2


Base.adjoint(g::EntanglementXXGate) = EntanglementXXGate(-g.θ[])
Base.adjoint(g::EntanglementYYGate) = EntanglementYYGate(-g.θ[])
Base.adjoint(g::EntanglementZZGate) = EntanglementZZGate(-g.θ[])


"""
General controlled gate
"""
struct ControlledGate{G} <: AbstractGate
    "target gate"
    U::G
    "number of control wires"
    M::Int
end

# wires
num_wires(g::ControlledGate)::Int = g.M + num_wires(g.U)
target_wires(g::ControlledGate) = num_wires(g.U)
control_wires(g::ControlledGate) = g.M

function matrix(g::ControlledGate{G}) where {G <:AbstractGate}
    N = num_wires(g)
    Umat = matrix(g.U)
    CU = Matrix{ComplexF64}(I, 2^N, 2^N)
    # Note: target qubit(s) corresponds to fastest varying index
    CU[end-size(Umat, 1)+1:end, end-size(Umat, 2)+1:end] = Umat
    return CU
end

function sparse_matrix(g::ControlledGate{G}) where {G <:AbstractGate}
    N = num_wires(g)
    Umat = sparse_matrix(g.U)
    CU = sparse(one(ComplexF64)*I, 2^N, 2^N)
    # Note: target qubit(s) corresponds to fastest varying index
    CU[end-size(Umat, 1)+1:end, end-size(Umat, 2)+1:end] = Umat
    return dropzeros!(CU)
end

LinearAlgebra.ishermitian(g::ControlledGate) = LinearAlgebra.ishermitian(g.U)

Base.adjoint(g::ControlledGate) = ControlledGate(Base.adjoint(g.U), g.M)

controlled_not() = ControlledGate(X, 1)


isunitary(m::AbstractMatrix) = (m * Base.adjoint(m) ≈ I)


"""
MatrixGate: general gate constructed from an unitary matrix
"""
struct MatrixGate <: AbstractGate
    matrix::SparseMatrixCSC{ComplexF64,Int}
    function MatrixGate(m::AbstractMatrix)
        @assert size(m, 1) == size(m, 2)
        isunitary(m) || error("Quantum gate must be unitary")
        return new(sparse(m))
    end
end

# wires
num_wires(g::MatrixGate)::Int = Int(log(2, size(g.matrix, 1)))

matrix(g::MatrixGate) = Matrix(g.matrix)

sparse_matrix(g::MatrixGate) = g.matrix

Base.adjoint(g::MatrixGate) = MatrixGate(Base.adjoint(g.matrix))

LinearAlgebra.ishermitian(g::MatrixGate) = LinearAlgebra.ishermitian(g.matrix)


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


# Define "addition" and "scaling" of parametric gates (required for gradient accumulation)

Base.:+(g1::RxGate, g2::RxGate) = RxGate(g1.θ[] + g2.θ[])
Base.:+(g1::RyGate, g2::RyGate) = RyGate(g1.θ[] + g2.θ[])
Base.:+(g1::RzGate, g2::RzGate) = RzGate(g1.θ[] + g2.θ[])

Base.:+(g1::RotationGate, g2::RotationGate) = RotationGate(g1.nθ + g2.nθ)

Base.:+(g1::PhaseShiftGate, g2::PhaseShiftGate) = PhaseShiftGate(g1.ϕ[] + g2.ϕ[])

Base.:+(g1::EntanglementXXGate, g2::EntanglementXXGate) = EntanglementXXGate(g1.θ[] + g2.θ[])
Base.:+(g1::EntanglementYYGate, g2::EntanglementYYGate) = EntanglementYYGate(g1.θ[] + g2.θ[])
Base.:+(g1::EntanglementZZGate, g2::EntanglementZZGate) = EntanglementZZGate(g1.θ[] + g2.θ[])

Base.:*(α::Real, g::RxGate) = RxGate(α * g.θ[])
Base.:*(α::Real, g::RyGate) = RyGate(α * g.θ[])
Base.:*(α::Real, g::RzGate) = RzGate(α * g.θ[])

Base.:*(α::Real, g::RotationGate) = RotationGate(α * g.nθ)

Base.:*(α::Real, g::PhaseShiftGate) = PhaseShiftGate(α * g.ϕ[])

Base.:*(α::Real, g::EntanglementXXGate) = EntanglementXXGate(α * g.θ[])
Base.:*(α::Real, g::EntanglementYYGate) = EntanglementYYGate(α * g.θ[])
Base.:*(α::Real, g::EntanglementZZGate) = EntanglementZZGate(α * g.θ[])
