# Models

```@meta
CurrentModule = Qaintessent
```

```@docs
qft_circuit(N::Integer)
```
```@docs
toffoli_circuit(cntrl::Tuple{<:Integer, <:Integer} , trg::Integer, N::Integer)
```

```@docs
vbe_adder_circuit(N::Integer)
```
#### Example
```jldoctest
using Qaintessent 

N = 3
M = N*3 + 1

adder = vbe_adder_circuit(N)
ψ = fill(0.0+0.0*im, 2^M)
a = 1
b = 3

index = b << N + a
ψ[index+1] = 1.0

ψ = apply(adder, ψ)
((findall(x->x==1, ψ)[1]-1)%(2^2N)) >> N

# output

4
```

```@docs
qcla_out_adder_circuit(N::Integer)
```
#### Example
```jldoctest
using Qaintessent 

N = 3
M = 3N + 1

adder = qcla_out_adder_circuit(N)
ψ = fill(0.0+0.0*im, 2^M)
a = 2
b = 3

index = b << N + a
ψ[index+1] = 1.0

ψ = apply(adder, ψ)
(findall(x->x==1, ψ)[1] - 1) >> 2N

# output

5
```

```@docs
qcla_inplace_adder_circuit(N::Integer)
```
#### Example
```jldoctest
using Qaintessent 

N = 3
M = 3N

adder = qcla_inplace_adder_circuit(N)
ψ = fill(0.0+0.0*im, 2^M)
a = 1
b = 2

index = b << N + a
ψ[index+1] = 1.0

ψ = apply(adder, ψ)
(findall(x->x==1, ψ)[1] - 1) >> N

# output

3
```

