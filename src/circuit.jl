

struct CircuitGate{M,N}
    "ordered wire indices which this gate acts on"
    iwire::NTuple{M, <:Integer}
    "actual gate"
    gate::AbstractGate{M}

    function CircuitGate{M,N}(iwire::NTuple{M, <:Integer}, gate::AbstractGate{M}) where {M,N}
        M ≥ 1 || error("Need at least one wire to act on.")
        M ≤ N || error("Number of gate wires cannot be larger than total number of wires.")
        length(unique(iwire)) == M || error("Wire indices must be unique.")
        minimum(iwire) ≥ 1 || error("Wire index cannot be smaller than 1.")
        maximum(iwire) ≤ N || error("Wire index cannot be larger than total number of wires.")
        new{M,N}(iwire, gate)
    end
end


function matrix(cg::CircuitGate{M,N}) where {M,N}
    # complementary wires
    iwcompl = setdiff(1:N, cg.iwire)
    @assert length(cg.iwire) + length(iwcompl) == N

    # TODO: support general "qudits"
    d = 2

    # TODO: handle sparse matrices efficiently
    gmat = matrix(cg.gate)
    @assert size(gmat) == (d^M, d^M)

    # Note: following the ordering convention of `kron` here, i.e.,
    # last qubit corresponds to fastest varying index
    strides = [d^(N-j) for j in 1:N]

    rowind = fill(0, d^(N+M))
    colind = fill(0, d^(N+M))
    values = fill(zero(eltype(gmat)), d^(N+M))
    count = 0
    for kw in (N-M > 0 ? reverse.(Iterators.product(fill(0:d-1, N-M)...)) : [()])
        for (i, iw) in enumerate(reverse.(Iterators.product(fill(0:d-1, M)...)))
            for (j, jw) in enumerate(reverse.(Iterators.product(fill(0:d-1, M)...)))
                # rearrange wire indices according to specification
                p = fill(0, N)
                for m in 1:(N-M); p[ iwcompl[m]] = kw[m]; end
                for m in 1:M;     p[cg.iwire[m]] = iw[m]; end
                q = fill(0, N)
                for m in 1:(N-M); q[ iwcompl[m]] = kw[m]; end
                for m in 1:M;     q[cg.iwire[m]] = jw[m]; end
                count += 1
                rowind[count] = dot(p, strides) + 1
                colind[count] = dot(q, strides) + 1
                values[count] = gmat[i, j]
            end
        end
    end
    @assert count == d^(N+M)

    return dropzeros!(sparse(rowind, colind, values, d^N, d^N))
end
