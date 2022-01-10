
"""
    AbstractGate

Abtract unitary quantum gate. `N` is the number of "wires" the gate acts on.
"""
abstract type AbstractGate end

"""
    Base.adjoint(::AbstractGate)

returns adjoint of AbstractGate objects

# Examples
```jldoctest    
julia> X ≈ X
true

julia> Y ≈ Z
false

julia> RxGate(0.3) ≈ RxGate(0.3)
true

julia> RxGate(0.3) ≈ RxGate(0.2)
false
```
"""
Base.adjoint(::AbstractGate) = error("Not implemented")

"""
    num_wires(::AbstractGate)

returns number of wires affected by AbstractGate object.

# Examples
```jldoctest
julia> num_wires(X)
1

julia> num_wires(SwapGate())
2
```
"""
num_wires(::AbstractGate) = error("Not implemented")

"""
    matrix(::AbstractGate)

returns dense matrix representation of [`AbstractGate`](@ref) objects

# Examples
```jldoctest    
julia> matrix(X)
2×2 Array{Complex{FloatQ},2}:
 0.0+0.0im  1.0+0.0im
 1.0+0.0im  0.0+0.0im
```
"""
matrix(::AbstractGate) = error("Not implemented")

"""
    sparse_matrix(::AbstractGate)

returns sparse matrix representation of [`AbstractGate`](@ref) objects

# Examples
```jldoctest    
julia> sparse_matrix(X)
2×2 SparseArrays.SparseMatrixCSC{Complex{FloatQ},Int64} with 2 stored entries:
  [2, 1]  =  1.0+0.0im
  [1, 2]  =  1.0+0.0im
```
"""
sparse_matrix(::AbstractGate) = error("Not implemented")

"""
    data(::AbstractGate)

returns variable data for specified gate

# Examples
```jldoctest    
julia> sparse_matrix(X)
```
"""
data(::AbstractGate) = nothing

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

const XMatrix = ComplexQ[0  1 ; 1  0]
const YMatrix = ComplexQ[0 -im; im 0]
const ZMatrix = ComplexQ[1  0 ; 0 -1]
const sXMatrix = sparse(XMatrix)
const sYMatrix = sparse(YMatrix)
const sZMatrix = sparse(ZMatrix)
matrix(::XGate) = XMatrix
matrix(::YGate) = YMatrix
matrix(::ZGate) = ZMatrix

sparse_matrix(::XGate) = sXMatrix
sparse_matrix(::YGate) = sYMatrix
sparse_matrix(::ZGate) = sZMatrix

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

``H = \\frac{1}{\\sqrt{2}} \\begin{pmatrix} 1 & 1 \\\\ 1 & 1 \\end{pmatrix}``
"""
struct HadamardGate <: AbstractGate end

const HMatrix = ComplexQ[1  1; 1 -1] ./ convert(ComplexQ, sqrt(2))
const sHMatrix = sparse(HMatrix) 

sparse_matrix(::HadamardGate) = sHMatrix
matrix(::HadamardGate) = HMatrix

LinearAlgebra.ishermitian(::HadamardGate) = true
# Hadamard gate is Hermitian
Base.adjoint(H::HadamardGate) = H

# corresponding instances
H = HadamardGate()
# wires
num_wires(::HadamardGate)::Int = 1


"""
    S gate

``S = \\begin{pmatrix} 1 & 0 \\\\ 0 & i \\end{pmatrix}``
"""
struct SGate <: AbstractGate end

"""
    T gate

``T = \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{\\frac{iπ}{4}} \\end{pmatrix}``
"""
struct TGate <: AbstractGate end

"""
    S^{†} gate

``S^{†} = \\begin{pmatrix} 1 & 0 \\\\ 0 & -i \\end{pmatrix}``
"""
struct SdagGate <: AbstractGate end

"""
    T^{†} gate

``T^{†} = \\begin{pmatrix} 1 & 0 \\\\ 0 & e^{-\\frac{iπ}{4}} \\end{pmatrix}``
"""
struct TdagGate <: AbstractGate end

const SGateMatrix    = ComplexQ[1 0; 0 im]
const SdagGateMatrix = ComplexQ[1 0; 0 -im]
const TGateMatrix    = ComplexQ[1 0; 0 Base.exp(im * π / 4)]
const TdagGateMatrix = ComplexQ[1 0; 0 Base.exp(-im * π / 4)]

const sSGateMatrix    = sparse(SGateMatrix)
const sSdagGateMatrix = sparse(SdagGateMatrix)
const sTGateMatrix    = sparse(TGateMatrix)
const sTdagGateMatrix = sparse(TdagGateMatrix)

matrix(::SGate)    = SGateMatrix
matrix(::SdagGate) = SdagGateMatrix
matrix(::TGate)    = TGateMatrix
matrix(::TdagGate) = TdagGateMatrix

sparse_matrix(::SGate)    = sSGateMatrix
sparse_matrix(::SdagGate) = sSdagGateMatrix
sparse_matrix(::TGate)    = sTGateMatrix
sparse_matrix(::TdagGate) = sTdagGateMatrix

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

``R_{x}(\\theta) = \\begin{pmatrix} \\cos(\\frac{\\theta}{2}) & -i\\sin(\\frac{\\theta}{2}) \\\\ -i\\sin(\\frac{\\theta}{2}) & \\cos(\\frac{\\theta}{2}) \\end{pmatrix}``
"""
struct RxGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{FloatQ}

    function RxGate(θ::Real)
        new([θ])
    end
end

const RxGateMatrix = ComplexQ[1 -im; -im 1]

function matrix(g::RxGate)
    m = zeros(ComplexQ, (2,2))
    c = cos(g.θ[] / 2.0)
    s = -im*sin(g.θ[] / 2.0)
    m[1,1] = c; m[1,2] = s;
    m[2,1] = s; m[2,2] = c;
    m
end

sparse_matrix(g::RxGate) = sparse(matrix(g))
data(g::RxGate) = g.θ[]
# data(::RxGate) = nothing
LinearAlgebra.ishermitian(g::RxGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::RxGate)::Int = 1


"""
    Rotation-Y gate

``R_{y}(\\theta) = \\begin{pmatrix} \\cos(\\frac{\\theta}{2}) & -\\sin(\\frac{\\theta}{2}) \\\\ \\sin(\\frac{\\theta}{2}) & \\cos(\\frac{\\theta}{2}) \\end{pmatrix}``
"""
struct RyGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{FloatQ}

    function RyGate(θ::Real)
        new([θ])
    end
end

function matrix(g::RyGate)
    m = zeros(ComplexQ, (2,2))
    c = cos(g.θ[] / 2.0)
    s = sin(g.θ[] / 2.0)
    m[1,1] = c; m[1,2] = -s;
    m[2,2] = c; m[2,1] = s
    m
end

sparse_matrix(g::RyGate) = sparse(matrix(g))
data(g::RyGate) = g.θ[]

LinearAlgebra.ishermitian(g::RyGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::RyGate)::Int = 1


"""
    Rotation-Z gate

``R_{z}(\\theta) = \\begin{pmatrix} e^{\\frac{-i\\theta}{2}} & 0 \\\\ 0 & e^{\\frac{i\\theta}{2}} \\end{pmatrix}``
"""
struct RzGate <: AbstractGate
    # use a reference type (array with 1 entry) for compatibility with Flux
    θ::Vector{FloatQ}
    function RzGate(θ::Real)
        new([θ])
    end
end

function matrix(g::RzGate)
    m = zeros(ComplexQ, (2,2))
    eθ = exp(im * g.θ[] / 2)
    m[1,1] = conj(eθ)
    m[2,2] = eθ
    m
end

sparse_matrix(g::RzGate) = sparse(matrix(g))
data(g::RzGate) = g.θ[]
LinearAlgebra.ishermitian(g::RzGate) = abs(sin(g.θ[] / 2)) < 4 * eps()


Base.adjoint(g::RxGate) = RxGate(-g.θ[])
Base.adjoint(g::RyGate) = RyGate(-g.θ[])
Base.adjoint(g::RzGate) = RzGate(-g.θ[])

# wires
num_wires(::RzGate)::Int = 1


"""
    General rotation operator gate: rotation by angle `θ` around unit vector `n`.

``R_{\\vec v}(\\theta) = \\cos(\\frac{\\theta}{2})I - i\\sin(\\frac{\\theta}{2})\\vec v \\sigma, \\\\ \\sigma = [X, Y, Z]``
"""
struct RotationGate <: AbstractGate
    nθ::Vector{FloatQ}

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
        return Matrix{ComplexQ}(I, 2, 2)
    end
    n = g.nθ / θ
    return cos(θ/2)*I - im*sin(θ/2)*pauli_vector(n...)
end

sparse_matrix(g::RotationGate) = sparse(matrix(g))
data(g::RotationGate) = to_gpu(matrix(g))
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
    ϕ::Vector{FloatQ}

    function PhaseShiftGate(ϕ::Real)
        new([ϕ])
    end
end

function matrix(g::PhaseShiftGate) 
    m = zeros(ComplexQ, (2,2))
    m[1,1] = 1; m[2,2] = Base.exp(im * g.ϕ[])
    m
end

sparse_matrix(g::PhaseShiftGate) = sparse(matrix(g))
data(g::PhaseShiftGate) = g.ϕ[]
LinearAlgebra.ishermitian(g::PhaseShiftGate) = abs(sin(g.ϕ[])) < 4 * eps()

Base.adjoint(g::PhaseShiftGate) = PhaseShiftGate(-g.ϕ[])

# wires
num_wires(::PhaseShiftGate)::Int = 1


"""
    Swap gate

``SWAP = \\begin{pmatrix} 1 & 0 & 0 & 0 \\\\ 0 & 0 & 1 & 0 \\\\ 0 & 1 & 0 & 0 \\\\ 0 & 0 & 0 & 1 \\end{pmatrix}``
"""
struct SwapGate <: AbstractGate end


const swapmatrix = ComplexQ[1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]
const sswapmatrix = sparse([1, 3, 2, 4], [1, 2, 3, 4], ComplexQ[1, 1, 1, 1])

matrix(::SwapGate,) = swapmatrix
sparse_matrix(::SwapGate,) = sswapmatrix

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
    θ::Vector{FloatQ}
    function EntanglementXXGate(θ::Real)
        new([θ])
    end
end

function matrix(g::EntanglementXXGate)
    m = zeros(ComplexQ, (4,4))
    c = cos(g.θ[] / 2)
    s = -im*sin(g.θ[] / 2)
    m[1,1] = c; m[1,4] = s;
    m[2,2] = c; m[2,3] = s;
    m[3,2] = s; m[3,3] = c;
    m[4,1] = s; m[4,4] = c;
    m
end

function sparse_matrix(g::EntanglementXXGate)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    sparse([1, 4, 2, 3, 2, 3, 1, 4], [1, 1, 2, 2, 3, 3, 4, 4], ComplexQ[c, -im*s, c, -im*s, -im*s, c, -im*s, c])
end

data(g::EntanglementXXGate) = g.θ[]

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
    θ::Vector{FloatQ}
    function EntanglementYYGate(θ::Real)
        new([θ])
    end
end

function matrix(g::EntanglementYYGate)
    m = zeros(ComplexQ, (4,4))
    c = cos(g.θ[] / 2)
    s = im*sin(g.θ[] / 2)
    m[1,1] = c; m[1,4] = s;
    m[2,2] = c; m[2,3] = -s;
    m[3,2] = -s; m[3,3] = c;
    m[4,1] = s; m[4,4] = c;
    m
end

function sparse_matrix(g::EntanglementYYGate)
    c = cos(g.θ[] / 2)
    s = sin(g.θ[] / 2)
    sparse([1, 4, 2, 3, 2, 3, 1, 4], [1, 1, 2, 2, 3, 3, 4, 4], ComplexQ[c, im*s, c, -im*s, -im*s, c, im*s, c])
end

data(g::EntanglementYYGate) = g.θ[]

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
    θ::Vector{FloatQ}
    function EntanglementZZGate(θ::Real)
        new([θ])
    end
end

function matrix(g::EntanglementZZGate)
    m = zeros(ComplexQ, (4,4))
    eθ = exp(im * g.θ[] / 2)
    m[1,1] = conj(eθ)
    m[2,2] = eθ
    m[3,3] = eθ
    m[4,4] = conj(eθ)
    m
end

function sparse_matrix(g::EntanglementZZGate)
    eθ = exp(im * g.θ[] / 2)
    sparse([1, 2, 3, 4], [1, 2, 3, 4], ComplexQ[conj(eθ), eθ, eθ, conj(eθ)])
end
data(g::EntanglementZZGate) = g.θ[]
LinearAlgebra.ishermitian(g::EntanglementZZGate) = abs(sin(g.θ[] / 2)) < 4 * eps()

# wires
num_wires(::EntanglementZZGate)::Int = 2


Base.adjoint(g::EntanglementXXGate) = EntanglementXXGate(-g.θ[])
Base.adjoint(g::EntanglementYYGate) = EntanglementYYGate(-g.θ[])
Base.adjoint(g::EntanglementZZGate) = EntanglementZZGate(-g.θ[])


matrix(::UniformScaling{Bool}) = ComplexQ[1 0; 0 1]
sparse_matrix(::UniformScaling{Bool}) = sparse((1.0+0.0im)*I, 2, 2)

"""
    kron(m::Union{UniformScaling{Bool},G}...) where {G<:AbstractGate}

returns sparse matrix representation of kronecker product of AbstractGate objects
"""
Base.kron(m::Union{UniformScaling{Bool},G}...) where {G<:AbstractGate} = Base.kron(sparse_matrix.(m)...)

"""
    ControlledGate{G}(U::G, M::Int) where {G<:AbstractGate}

returns ControlledGate object. Controls [`AbstractGate`](@ref) of type `G` with `M` control wires
"""
struct ControlledGate{G} <: AbstractGate
    "target gate"
    U::G
    "number of control wires"
    M::Int
end

# wires
num_wires(g::ControlledGate)::Int = g.M + num_wires(g.U)

"""
    target_wires(g::ControlledGate)

return target wires for a [`ControlledGate`](@ref)

# Examples
```jldocttest
# controlled gate with 1 target wire
julia> cnot = ControlledGate(X, 1);
julia> target_wires(cnot)
1
```
"""
target_wires(g::ControlledGate) = num_wires(g.U)
target_wires(g::AbstractGate) = num_wires(g)

"""
    control_wires(g::ControlledGate)

return control wires for a [`ControlledGate`](@ref)

# Examples
```jldocttest
# controlled gate with 2 control wires
julia> ccnot = ControlledGate(X, 2);
julia> target_wires(ccnot)
2
```
"""
control_wires(g::ControlledGate) = g.M
control_wires(g::AbstractGate) = 0

function matrix(g::ControlledGate{G}) where {G <:AbstractGate}
    N = num_wires(g)
    Umat = matrix(g.U)
    CU = Matrix{ComplexQ}(I, 2^N, 2^N)
    # Note: target qubit(s) corresponds to fastest varying index
    CU[end-size(Umat, 1)+1:end, end-size(Umat, 2)+1:end] = Umat
    return CU
end

function sparse_matrix(g::ControlledGate{G}) where {G <:AbstractGate}
    N = num_wires(g)
    Umat = sparse_matrix(g.U)
    CU = sparse(one(ComplexQ)*I, 2^N, 2^N)
    # Note: target qubit(s) corresponds to fastest varying index
    CU[end-size(Umat, 1)+1:end, end-size(Umat, 2)+1:end] = Umat
    return dropzeros!(CU)
end
data(g::ControlledGate) = data(g.U)
LinearAlgebra.ishermitian(g::ControlledGate) = LinearAlgebra.ishermitian(g.U)

Base.adjoint(g::ControlledGate) = ControlledGate(Base.adjoint(g.U), g.M)

controlled_not() = ControlledGate(X, 1)


isunitary(m::AbstractMatrix) = (m * Base.adjoint(m) ≈ I)


"""
    MatrixGate

general gate constructed from an unitary matrix
"""
struct MatrixGate <: AbstractGate
    matrix::SparseMatrixCSC{ComplexQ,Int}
    function MatrixGate(m::AbstractMatrix)
        @assert size(m, 1) == size(m, 2)
        isunitary(m) || error("Quantum gate must be unitary")
        return new(sparse(m))
    end
end

# wires
num_wires(g::MatrixGate)::Int = Int(log(2, size(g.matrix, 1)))

matrix(g::MatrixGate) = Matrix{ComplexQ}(g.matrix)

sparse_matrix(g::MatrixGate) = g.matrix

data(g::MatrixGate) = to_gpu(g.matrix)

Base.adjoint(g::MatrixGate) = MatrixGate(Base.adjoint(g.matrix))

LinearAlgebra.ishermitian(g::MatrixGate) = LinearAlgebra.ishermitian(g.matrix)

"""
    Base.isapprox(g1::G, g2::G) where {G <: AbstractGate} 

returns true if input AbstractGates are of same type with same paremeters

# Examples
```jldoctest
julia> X ≈ X
true

julia> X ≈ Y
false

julia> RxGate(0.1) ≈ RxGate(0.2)
false
```
"""
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

