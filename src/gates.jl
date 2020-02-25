
"""
    AbstractGate{N}

Abtract unitary quantum gate. 'N' is the number of "wires" the gate acts on.
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


struct HadamardGate <: AbstractGate{1} end

matrix(::HadamardGate) = [1 1; 1 -1] / sqrt(2)

# Hadamard gate is Hermitian
Base.adjoint(H::HadamardGate) = H


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
