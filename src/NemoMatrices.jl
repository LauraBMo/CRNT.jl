
export Diagonal, Jacobian

"""

    matrix(R, V::AbstractVector{T}) where {T}

Given a Julia vector `V` of entries, construct the corresponding AbstractAlgebra.jl one-column matrix over the given ring `R`, assuming all the entries can be coerced into `R`.
"""
function Nemo.matrix(R, V::AbstractVector{T}) where {T}
    return Nemo.matrix(R, reshape(V, (:,1)))
end

"""

    Diagonal(R, V::AbstractVector{T}) where {T}

Given a Julia vector `V` of entries, construct the corresponding AbstractAlgebra.jl diagonal matrix over the given ring `R`, assuming all the entries can be coerced into `R`.
"""
function Diagonal(R, V::AbstractVector{T}) where {T}
    n = size(V,1)
    D = Nemo.zero_matrix(R,n,n)
    for i in 1:n D[i,i] = V[i] end
    return D
end

"""

    Jacobian(R, M)
    Jacobian(R, M, xs=1:length(Nemo.gens(R)))

Returns the Jacobian matrix of a one-column matrix of polynomials `M` with respect to the generators of `R` indexed by `xs`. When `xs` is omitted all the generators of `R` are used.

# Examples
```jldoctest; setup = :(using CRNT)
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> M = matrix(R, [vars[1]*vars[2]*vars[6]-vars[8]; vars[3]*vars[9]+2*vars[7]])
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
function Jacobian(R, M, xs=1:length(Nemo.gens(R)))
    r, c = size(M)
    @assert c == 1
    J = Nemo.zero_matrix(R, r, length(xs))
    for i in 1:r
        for (j, k) in enumerate(xs)
            J[i,j] = Nemo.derivative(M[i,1], k)
        end
    end
    return J
end
