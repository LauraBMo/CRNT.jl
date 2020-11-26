
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

function evenodd_coeffs(p)
    l = length(p); ## degree +1
    ## For the zero polynomial, 'odd' would be a 0-element array
    if l < 1
        even = [0]
        odd  = [0]
    else
        even = [coeff(p, i) for i in 0:2:l]
        odd  = [coeff(p, i) for i in 1:2:l]
    end
    return even, odd
end

function evenodd_polys(p)
    even, odd = evenodd_coeffs(p)
    return parent(p)(evalxsqrt(even)), parent(p)(evalxsqrt(odd))
end

function evalxsqrt(poly)
    n = size(poly, 1)
    R = parent(poly[1])
    out = fill(R(0), 2 * n - 1)
    out[1] = poly[1]
    for i in 2:n
        out[2 * i - 1] = poly[i]
    end
    return out
end

function HurwitzMatrix1(R::AbstractAlgebra.Ring, p)
    sylvester_matrix(evenodd_polys(p)...)
end


# ## k's are parameters
# R, (k1, k2) = PolynomialRing(Nemo.QQ, ["k1", "k2"]);
# RQ = FractionField(R)

# ## The ring for the polynomial system
# S, (x, y, z, w, t) = PolynomialRing(RQ, ["x", "y", "z", "w", "t"])
# SQ = FractionField(S)

# ## The ring for the characteristic polynomial.
# T, l = PolynomialRing(SQ, "A")

# P = T([x * y * k1, k1 * k2^2 + 2, k1, //(k1, k2 + k1 * k2) * x])
# Q = T([x, y, z, w, t])
# p, q = evenodd_polys(Q)

# sylvester_matrix(evenodd_polys(p)...)
# -x * resultant_subresultant(p, q) == det(HurwitzMatrix(S, Q))

# # add_row(a::MatrixElem, s::RingElement, i::Int, j::Int, cols = 1:ncols(a))

# R, (k1, k2, k3) = PolynomialRing(Nemo.QQ, ["k1", "k2", "k3"]);
# RQ = FractionField(R)

# ## The ring for the polynomial system
# S, (x, y, z) = PolynomialRing(RQ, ["x", "y", "z"])

# ## The ring for the characteristic polynomial.
# T, l = PolynomialRing(S, "A")

# chp = T([k1 * k2 * y + k2 * k3 * z + k1 * k3, k2 * y + k2 * z + k1 + k3, S(1)])

# p, q = evenodd_polys(chp)
# resultant_subresultant(p,q)
# M = sylvester_matrix(p, q)

# det(M) == resultant_subresultant(p, q)

# det(M) == det(HurwitzMatrix(S, chp))

function ISRO!(M, j)
    c = ncols(M)
    p = M[j,j]
    for i in (j + 1):c
        if !(iszero(M[i,j]))
            add_row!(M, M[i,j] // p, j, i, c)
        end
    end
end

function ISRO!(M)
    return j -> ISRO!(M, j)
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
        # signs = unique(sign.(coeffs(D)))
        signs = unique(sign.(vec(reshape(hcat, coeffs.(coeffs(D))))))
        write(io, "$(signs)\n")
        if -1 in signs
            H = [H; D]
        end
    end
    write(io, "==============================================\n")
    write(io, "==============================================\n")
    return H
end
