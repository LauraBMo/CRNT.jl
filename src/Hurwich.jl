
struct HurwitzDeterminants{T}
    M::MatElem{T}
end

Base.convert(::Type{AbstractMatrix},M::HurwitzDeterminants) = Array(M.M)
Base.show(io::IO, M::HurwitzDeterminants) = print(io, M.M)
Base.show(io::IO, ::MIME"text/plain", z::HurwitzDeterminants{T}) where {T} =
    print(io, "HurwitzDeterminants{$T} iterator:\n", z)

function Base.iterate(M::HurwitzDeterminants, state=1)
    if state > size(M.M, 1)
        return nothing
    else
        return Nemo.det(M.M[1:state,1:state]), state + 1
    end
end

# https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-iteration-1
Base.length(M::HurwitzDeterminants) = size(M.M, 1)
Base.eltype(::Type{HurwitzDeterminants{T}}) where {T} = T

# Input: polynomial p in one variable
function HurwitzMatrix(S, p)
    n = length(p) - 1
    M = zero_matrix(S, n, n)
    for i in 1:n
        for j in 1:n
            if (n - 2 * i + j) in 0:n
                M[i,j] = coeff(p, n - 2 * i + j)
            end
        end
    end
    return M
end

function HurwitzDeterminants(S, p)
    return HurwitzDeterminants(HurwitzMatrix(S, p))
end


# mindeg = size(N,1)-rank(matrix(S,N))
function StabilityMatrix(S, J, mindeg)
    P, x = PolynomialRing(S, "x")
    p = divexact(charpoly(P, J), x^mindeg)
    iter = HurwitzDeterminants(S, p)
    print("$(length(iter)) expressions will be studied")
    for (i, D) in enumerate(iter)
        print("$i th determinant")
        print("$(unique(sign.(coeffs(D))))")
    end
end
