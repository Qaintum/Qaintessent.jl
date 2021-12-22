using Base
using Base: invpermute!!
using Qaintessent
using CUDA

"""Permutes perm vector in Statevector object, shifting qubits in "wires" to the fastest running qubit"""
function _orderbits(a, b, c, iwire, N)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        result = 0
        r = length(iwire)
        w = 0
        for n in 1:N
            bit = ((i-1) >> (n-1)) & 1
            if n in iwire
                result = result | (bit << w)
                w = w + 1
            else
                result = result | (bit << r)
                r = r + 1
            end
        end
        @inbounds c[result+1] = a[i]
        @inbounds b[i] = result+1
    end
    return
end

"""Permutes perm vector in Statevector object, shifting qubits in "wires" to the fastest running qubit"""
function _permutebits(a, b, c)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        @inbounds c[i] = a[b[i]]
    end
    return
end

for g = [:XGate, :YGate, :ZGate, :SGate, :SdagGate, :TGate, :TdagGate, :RxGate, :RyGate, :RzGate, :PhaseShiftGate, :HadamardGate, :RotationGate, :MatrixGate]
    func_name = Symbol("_apply", g)
    eval(quote
        """Tailored gpu apply for g"""
        function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::CircuitGate{1,$g})
            iwire = cg.iwire[1]
            bin = 1 << (iwire-1)
            CUDA.@sync @cuda threads=1024 blocks=length(a)÷1024+1 $func_name(a, b, c, iwire, bin, data(cg))
        end
    end)

    eval(quote
        """Tailored gpu apply for g"""
        function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::MeasurementOperator{1,$g})
            iwire = cg.iwire[1]
            bin = 1 << (iwire-1)
            CUDA.@sync @cuda threads=1024 blocks=length(a)÷1024+1 $func_name(a, b, c, iwire, bin, data(cg))
        end
    end)
end

for g = [:SwapGate]
    func_name = Symbol("_apply", g)
    eval(quote
        """Tailored gpu apply for g"""
        function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::CircuitGate{2,$g})
            bin1 = 1 << (cg.iwire[1]-1)
            bin2 = 1 << (cg.iwire[2]-1)
            CUDA.@sync @cuda threads=1024 blocks=length(a)÷1024+1 $func_name(a, b, c, cg.iwire[1], cg.iwire[2], bin1, bin2, data(cg))
        end
    end)

    eval(quote
    """Tailored gpu apply for g"""
    function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::MeasurementOperator{2,$g})
        iwire = cg.iwire[1]
        bin = 1 << (iwire-1)
        CUDA.@sync @cuda threads=1024 blocks=length(a)÷1024+1 $func_name(a, b, c, iwire, bin, data(cg))
    end
end)

end

for g = [:XGate, :YGate, :ZGate, :SGate, :SdagGate, :TGate, :TdagGate, :RxGate, :RyGate, :RzGate, :PhaseShiftGate, :HadamardGate, :RotationGate]
    func_name = Symbol("_apply", g)
    eval(quote
        """Tailored gpu apply for g"""
        function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::CircuitGate{M,ControlledGate{$g}}) where {M}
            icontrol = cg.iwire[2:M]
            iwire = cg.iwire[1]
            mask = 0
            for i in icontrol
                mask += (1 << (i-1))
                if i < cg.iwire[1]
                    iwire -= 1
                end
            end
            d = a[(b.-1) .& mask .== mask]
            e = b[1:length(d)]
            f = c[(b.-1) .& mask .== mask]
            bin = 1 << (iwire-1)

            CUDA.@sync @cuda threads=1024 blocks=length(d)÷1024+1 $func_name(d, e, f, iwire, bin, data(cg))
            c .= a
            c[(b.-1) .& mask .== mask] .= f
            return
        end
    end)
end

for g = [:SwapGate]
    func_name = Symbol("_apply", g)
    eval(quote
        """Tailored gpu apply for g"""
        function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::CircuitGate{M,ControlledGate{$g}}) where {M}
            icontrol = cg.iwire[3:M]
            iwire1 = cg.iwire[1]
            iwire2 = cg.iwire[2]
            mask = 0
            for i in icontrol
                mask += (1 << (i-1))
                if i < cg.iwire[1]
                    iwire1 -= 1
                end
                if i < cg.iwire[2]
                    iwire2 -= 1 
                end
            end
            d = a[(b.-1) .& mask .== mask]
            e = b[1:length(d)]
            f = c[1:length(d)]

            bin1 = 1 << (iwire1 - 1)
            bin2 = 1 << (iwire2 - 1)
            CUDA.@sync @cuda threads=1024 blocks=length(d)÷1024+1 $func_name(d, e, f, iwire1, iwire2, bin1, bin2, data(cg))
            c .= a
            c[(b.-1) .& mask .== mask] .= f
            return
        end
    end)
end

"""Tailored gpu apply for general controlled CircuitGate"""
function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::CircuitGate{M,ControlledGate{G}}) where {M,G<:AbstractGate}
    U = CuArray(reshape(transpose(sparse_matrix(cg.gate.U)), :))
    T = Qaintessent.target_wires(cg.gate)
    C = Qaintessent.control_wires(cg.gate)
    N = 

    iwire = convert.(Int32, collect(cg.iwire[1:T]))
    icontrol = cg.iwire[T+1:M]
    mask = 0
    for i in icontrol
        mask += (1 << (i-1))
    end

    for (j,k) in enumerate(iwire)
        iwire[j] -= sum(icontrol .< (iwire[j],))
    end
    d = a[(b.-1) .& mask .== mask]
    e = b[1:length(d)]
    f = c[(b.-1) .& mask .== mask]
    iwire = CuArray(iwire)
    T = convert(Int32, T)
    N = convert(Int32, log2(length(a))-C)
    
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _orderbits(d, e, f, iwire, N) end
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _applyGate(d, e, f, T, U) end
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _permutebits(d, e, f) end
    c .= a
    c[(b.-1) .& mask .== mask] .= f
    return
end

"""Tailored gpu apply for general CircuitGate"""
function _apply!(a::CuArray, b::CuArray, c::CuArray, cg::CircuitGate{M,G}) where {M,G<:AbstractGate}
    U = CuArray(reshape(transpose(sparse_matrix(cg.gate)), :))
    iwire = CuArray(convert.(Int32, collect(cg.iwire)))
    T = convert(Int32, M)
    N = convert(Int32, log2(length(a)))
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _orderbits(a, b, c, iwire, N) end
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _applyGate(a, b, c, T, U) end
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _permutebits(a, b, c) end
    return
end

"""Tailored gpu apply for general MeasurementOperator"""
function _apply!(a::CuArray, b::CuArray, c::CuArray, mea::MeasurementOperator{M,G}) where {M,G<:Union{AbstractGate,AbstractMatrix}}
    U = CuArray(reshape(transpose(sparse_matrix(mea.operator)), :))
    iwire = CuArray(convert.(Int32, collect(mea.iwire)))
    num_wires = convert(Int32, M)
    system_size = convert(Int32, log2(length(a)))
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _orderbits(a, b, c, iwire, system_size) end
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _applyGate(a, b, c, num_wires, U) end
    CUDA.@sync begin @cuda threads=1024 blocks=length(a)÷1024+1 _permutebits(a, b, c) end
    return
end


"""General apply for general gates"""
function _applyGate(a, b, c, wires, U)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        a[i] = 0.0+0.0im
        bin = 1 << wires
        step = (i-1) % bin
        start = ((i-1) >> wires) << wires
        for j in 0:(bin-1)
            @inbounds a[i] += U[step*bin+j+1]*c[start+j+1]
        end
        return
    end
end


"""Tailored apply for XGate"""
function _applyXGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        @inbounds c[i] = a[(b[i] - 1) ⊻ bin + 1]
    end
    return
end


"""Tailored apply for YGate"""
function _applyYGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    pow = ((i-1) >> (wire-1)) & 1
    if i <= length(a)
        @inbounds c[i] = -im*a[(b[i] - 1) ⊻ bin + 1] * (-1)^pow
    end
    return
end


"""Tailored apply for ZGate"""
function _applyZGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    pow = ((i-1) >> (wire-1)) & 1
    if i <= length(a)
        @inbounds c[i] = a[i] * (-1)^pow
    end
    return
end


"""Tailored apply for HadamardGate"""
function _applyHadamardGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    pow = ((i-1) >> (wire-1)) & 1
    if i <= length(a)
        @inbounds c[i] = ( (-1)^pow*a[i] + a[(b[i] - 1) ⊻ bin + 1])/sqrt(2)
    end
    return
end


"""Tailored apply for SGate"""
function _applySGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x

    if i <= length(a)
        if (i-1) & bin == bin
            @inbounds c[i] = im*a[i]
        else 
            @inbounds c[i] = a[i]
        end
    end
    return
end


"""Tailored apply for SdagGate"""
function _applySdagGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x

    if i <= length(a)
        if (i-1) & bin == bin
            @inbounds c[i] = -im*a[i]
        else 
            @inbounds c[i] = a[i]
        end
    end
    return
end


"""Tailored apply for TGate"""
function _applyTGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    eπ4 = Base.exp(im*π/4)
    if i <= length(a)
        if (i-1) & bin == bin
            @inbounds c[i] = eπ4*a[i]
        else 
            @inbounds c[i] = a[i]
        end
    end
    return
end


"""Tailored apply for TdagGate"""
function _applyTdagGate(a, b, c, wire, bin, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    eπ4 = Base.exp(-im*π/4)
    if i <= length(a)
        if (i-1) & bin == bin
            @inbounds c[i] = eπ4*a[i]
        else 
            @inbounds c[i] = a[i]
        end
    end
    return
end


"""Tailored apply for RxGate"""
function _applyRxGate(a, b, c, wire, bin, θ)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    cθ = cos(θ/2)
    sθ = -im*sin(θ/2)
    if i <= length(a)
        j = (b[i] - 1) ⊻ bin + 1
        if (i-1) & bin == bin
            @inbounds c[i] = cθ*a[i] + sθ*a[j]
        else
            @inbounds c[i] = sθ*a[j] + cθ*a[i]
        end
        return
    end
    return
end


"""Tailored apply for RyGate"""
function _applyRyGate(a, b, c, wire, bin, θ)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    cθ = cos(θ/2)
    sθ = sin(θ/2)
    if i <= length(a)
        j = (b[i] - 1) ⊻ bin + 1
        if (i-1) & bin != bin
            @inbounds c[i] = cθ*a[i] - sθ*a[j]
        else
            @inbounds c[i] = sθ*a[j] + cθ*a[i]
        end
        return
    end
    return
end


"""Tailored apply for RzGate"""
function _applyRzGate(a, b, c, wire, bin, θ)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    neπ2 = Base.exp(-im*θ/2)
    eπ2 = Base.exp(im*θ/2)
    if i <= length(a)
        if (i-1) & bin != bin
            @inbounds c[i] = neπ2*a[i]
        else 
            @inbounds c[i] = eπ2*a[i]
        end
    end
    return
end


"""Tailored apply for PhaseShiftGate"""
function _applyPhaseShiftGate(a, b, c, wire, bin, ϕ)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    eϕ = Base.exp(im*ϕ)
    if i <= length(a)
        if (i-1) & bin == bin
            @inbounds c[i] = eϕ*a[i]
        else
            @inbounds c[i] = a[i]
        end
    end
    return
end

"""Tailored apply for RotationGate"""
function _applyRotationGate(a, b, c, wire, bin, U)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        if (i-1) & bin != bin
            @inbounds c[i] = U[1,1]*a[i] + U[1,2]*a[(b[i] - 1) ⊻ bin + 1]
        else
            @inbounds c[i] = U[2,1]*a[(b[i] - 1) ⊻ bin + 1] + U[2,2]*a[i]
        end
        return
    end
    return
end


"""Tailored apply for a general single qubit gate"""
function _applyMatrixGate(a, b, c, wire, bin, U)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        if (i-1) & bin != bin
            @inbounds c[i] = U[1,1]*a[i] + U[1,2]*a[(b[i] - 1) ⊻ bin + 1]
        else
            @inbounds c[i] = U[2,1]*a[(b[i] - 1) ⊻ bin + 1] + U[2,2]*a[i]
        end
        return
    end
    return
end


"""Tailored apply for SwapGate"""
function _applySwapGate(a, b, c, wire1, wire2, bin1, bin2, data=nothing)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x - 1
    bit1 = (i >> (wire1-1)) & 1
    bit2 = (i >> (wire2-1)) & 1
    x = (bit1 ⊻ bit2)
    x = (x << (wire1-1)) | (x << (wire2-1))
    result = i ⊻ x
    if i < length(a)
        @inbounds c[i+1] = a[result+1]
    end
    return
end


"""
    apply(ψ::Statevector, cg::CircuitGate{M,G}) where {M,G}

Apply a [`CircuitGate`](@ref) to a quantum state vector `ψ`.
# Examples

```jldoctest
julia> cg = circuit_gate(1, HadamardGate());
julia> ψ = Statevector(ComplexF64[1 0]);
julia> apply!(ψ, cg)
julia> ψ.state
2-element Array{Complex{Float64},1}:
 0.7071067811865475 + 0.0im
 0.7071067811865475 + 0.0im
```
"""
function apply!(ψ::Statevector, cg::CircuitGate{M,G}) where {M,G<:AbstractGate}
    req_wires(cg) <= ψ.N || error("CircuitGate requires a minimum of $(req_wires(cg)) qubits, input vector `ψ` has $ψ.N qubits")
    vector = CuArray(ψ.state)
    perm = CuArray(ψ.perm)
    scratch = CuArray(ψ.vec)
    _apply!(vector, perm, scratch, cg)
    ψ.state .= Array(scratch)
    return
end


"""
    apply(ψ::Statevector, cgs::Vector{<:CircuitGate})

Apply a sequence of [`CircuitGate`](@ref)(s) to a quantum state vector `ψ`.

# Examples

```jldoctest
julia> cgs = [circuit_gate(1, HadamardGate()),
                circuit_gate(1, X),
                circuit_gate(1, Y)];
julia> ψ = Statevector(ComplexF64[1 0]);
julia> apply!(ψ, cgs)
julia> ψ.state
2-element Array{Complex{Float64},1}:
 0.0 - 0.7071067811865475im
 0.0 + 0.7071067811865475im
```
"""
function apply!(ψ::Statevector{M}, cgs::Vector{<:CircuitGate}) where {M<:Complex}
    maximum(req_wires.(cgs)) <= ψ.N || error("CircuitGates require a minimum of $req qubits, input vector `ψ` has $N qubits")
    vector = CuArray(ψ.state)
    perm = CuArray(ψ.perm)
    scratch = CuArray(ψ.vec)
    for cg in cgs 
        _apply!(vector, perm, scratch, cg)
        tmp = vector
        vector = scratch
        scratch = tmp
    end
    ψ.state .= Array(vector)
    return
end


"""
    apply(ψ::Statevector, m::Moment)

returns state vector of `N` qubits after applying a `Moment{N}` object to a quantum state vector of `N` qubits `ψ`
"""
function apply!(ψ::Statevector, m::Moment) 
    Nmoment = req_wires(m)
    Nmoment <= ψ.N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    vector = CuArray(ψ.state)
    perm = CuArray(ψ.perm)
    scratch = CuArray(ψ.vec)
    for gate in m
        _apply!(vector, perm, scratch, m)
        tmp = vector
        vector = scratch
        scratch = tmp
    end
    ψ.state .= Array(vector)
    return
end


function _apply!(vector::CuArray, perm::CuArray, scratch::CuArray, m::Moment)
    for gate in m
        _apply!(vector, perm, scratch, gate)
        tmp = vector
        vector = scratch
        scratch = tmp
    end
    return vector, scratch
end


function apply!(ψ::Statevector, m::Vector{Moment})
    length(m) != 0 || error("Vector of length 0 cannot be applied")
    Nmoment = maximum(req_wires.(m))
    Nmoment <= ψ.N || error("Moment affecting $Nmoment qubits applied to $N qubits")
    vector = CuArray(ψ.state)
    perm = CuArray(ψ.perm)
    scratch = CuArray(ψ.vec)
    for moment in m
        vector, scratch = _apply!(vector, perm, scratch, moment)
    end
    ψ.state .= Array(vector)
    return
end


"""
    expectation_value!(vector, scratch, meas)

calculates expectation value based on output and measurement
"""
function _expectation_value(a, b, output)
    i = (blockIdx().x-1) * blockDim().x + threadIdx().x
    if i <= length(a)
        @inbounds output[i] = conj(a[i]) * b[i]
    end
    return
end


"""
    apply(ψ::Statevector, c::Circuit{N}) where {N}

returns list of expectation values from measurement operators in `c.meas` after applying circuit gates in `c.cgc` on state vector of `N` qubits `ψ`
"""
function apply!(ψ::Statevector{M}, c::Circuit{N}) where {M<:Complex,N}
    N == ψ.N || error("Size of vector `ψ` must match Circuit size of $(2^N)")
    length(c) != 0 || error("Circuit does not contain any gates")
    length(c.meas) != 0 || error("Circuit does not contain any measurement operators")
    vector = CuArray(ψ.state)
    perm = CuArray(ψ.perm)
    scratch = CuArray(ψ.vec)
    stored_state = similar(vector)
    agg = similar(vector)
    output = zeros(FloatQ, length(c.meas))

    for moment in c.moments
        vector, scratch = _apply!(vector, perm, scratch, moment)
    end
    ψ.state .= Array(vector)
    ψ.vec .= Array(scratch)
    stored_state .= vector

    for i in 1:length(c.meas)
        _apply!(vector, perm, scratch, c.meas[i])
        vector .= stored_state
        CUDA.@sync @cuda threads=1024 blocks=length(vector)÷1024+1 _expectation_value(vector, scratch, agg)
        CUDA.@sync output[i] = real(reduce(+, agg))
    end
    return Array(output)
end