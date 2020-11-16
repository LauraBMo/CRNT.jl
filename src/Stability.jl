
export HurwitzMatrix,
    HurwitzDeterminants,
    StabilityMatrix


struct HurwitzDeterminants{T}
    M::MatElem{T}
end

Base.convert(::Type{AbstractMatrix},M::HurwitzDeterminants) = Array(M.M)
Base.show(io::IO, M::HurwitzDeterminants) = print(io, M.M)
Base.show(io::IO, ::MIME"text/plain", z::HurwitzDeterminants{T}) where {T} =
    print(io, "HurwitzDeterminants{$T} iterator:\n", z)

# function Base.iterate(M::HurwitzDeterminants)
#     state = 1
#     return Base.iterate(M, state)
# end

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
function HurwitzMatrix(R::AbstractAlgebra.Ring, p)
    n = length(p) - 1
    M = zero_matrix(R, n, n)
    for i in 1:n
        for j in 1:n
            if (n - 2 * i + j) in 0:n
                M[i,j] = coeff(p, n - 2 * i + j)
            end
        end
    end
    return M
end

function HurwitzDeterminants(R::AbstractAlgebra.Ring, p)
    return HurwitzDeterminants(HurwitzMatrix(R, p))
end

# mindeg = size(N,1)-rank(matrix(R,N))
function StabilityMatrix(R::AbstractAlgebra.Ring, J, mindeg::Integer=size(J, 1) - rank(J))
    P, x = PolynomialRing(R, "x")
    p = divexact(charpoly(P, J), x^mindeg)
    H = []
    iter = HurwitzDeterminants(R, p)
    print("==============================================\n")
    print("==============================================\n")
    print("$(length(iter)) expressions will be studied\n")
    for (i, D) in enumerate(iter)
        print("\n-- $i th determinant --\n")
        signs = unique(sign.(coeffs(D)))
        print("$(signs)\n")
        if -1 in signs
            H = [H; D]
        end
    end
    print("==============================================\n")
    print("==============================================\n")
    return H
end

function StabilityMatrix(io, R::AbstractAlgebra.Ring, J, mindeg::Integer=size(J, 1) - rank(J))
    P, x = PolynomialRing(R, "x")
    p = divexact(charpoly(P, J), x^mindeg)
    H = []
    iter = HurwitzDeterminants(R, p)
    write(io, "==============================================\n")
    write(io, "==============================================\n")
    write(io, "$(length(iter)) expressions will be studied\n")
    for (i, D) in enumerate(iter)
        write(io, "\n-- $i th determinant --\n")
        signs = unique(sign.(coeffs(D)))
        write(io, "$(signs)\n")
        if -1 in signs
            H = [H; D]
        end
    end
    write(io, "==============================================\n")
    write(io, "==============================================\n")
    return H
end
