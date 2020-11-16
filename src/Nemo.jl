
export dissect,
    Diagonal,
    partialdervativeof,
    Jacobian

@doc raw"""
    dissect(p::MPolyElem)

Iterator for the coefficients and exponent vectors of the given polynomial.

"""
function dissect(p::MPolyElem)
    return zip(coeffs(p), exponent_vectors(p))
end

function Nemo.nullspace_right_rational(N::AbstractArray{T}) where {T<:Integer}
    Nnemo = Nemo.matrix(Nemo.FlintZZ, N)
    r, U = Nemo.nullspace_right_rational(Nnemo)
    return convert_to_array(U[:,1:r],T(1))
end

"""

    matrix(R, V::AbstractVector{T}) where {T}

Given a Julia vector `V` of entries, construct the corresponding AbstractAlgebra.jl one-column matrix over the given ring `R`, assuming all the entries can be coerced into `R`.
"""
function Nemo.matrix(R::AbstractAlgebra.Ring, V::AbstractVector)
    return Nemo.matrix(R, reshape(V, (:, 1)))
end


@doc raw"""

    Diagonal(R, V::AbstractVector{T}) where {T}

Given a Julia vector `V` of entries, construct the corresponding AbstractAlgebra.jl diagonal matrix over the given ring `R`, assuming all the entries can be coerced into `R`.
"""
function Diagonal(R::AbstractAlgebra.Ring, V::AbstractVector{T}) where {T}
    n = size(V, 1)
    D = Nemo.zero_matrix(R, n, n)
    for i in 1:n D[i,i] = V[i] end
    return D
end

function partialdervativeof(i::Integer)
    return p -> Nemo.derivative(p, i)
end

"""

    Jacobian(R, M)
    Jacobian(R, M, vars=1:length(gens(R)))

Returns the Jacobian matrix of a one-column matrix of polynomials `M` with respect to the generators of `R` indexed by `vars`. When `vars` is omitted all the generators of `R` are used.

# Examples
```jldoctest; setup = :(using CRNT, Nemo)
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> M = [vars[1]*vars[2]*vars[6]-vars[8]; vars[3]*vars[9]+2*vars[7]]
[k1*k2*x1-x3]
[ k3*x4+2*x2]

julia> Jacobian(R, M, 6:9)
[k1*k2  0  -1   0]
[    0  2   0  k3]

julia> Jacobian(R,M)
[k2*x1  k1*x1   0  0  0  k1*k2  0  -1   0]
[    0      0  x4  0  0      0  2   0  k3]
```
"""
function Jacobian(R::MPolyRing, M::AbstractVector, vars=1:length(gens(R)))
    rows = size(M, 1)
    cols = length(vars)
    J = fill(zero(R), rows, cols)
    for j in 1:cols
        J[:,j] = partialdervativeof(vars[j]).(M)
    end
    return J
end
