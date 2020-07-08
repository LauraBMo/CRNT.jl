module CRNT

# Write your package code here.

export NonnegativeNullspaceCone

# We need Nemo for the kernel of N1
using Nemo # matrix, FlintIntegerRing, nullspace
import AbstractAlgebra: Ring
import Polymake # polytope.cone,intersection
import LinearAlgebra: I, dot
import PolynomialRoots: roots

## Nemo: matrix, FlintZZ, nullspace
## Polymake: polytope:cone,intersection
## LinearAlgebra: I
"""

    NonnegativeNullspaceCone(N::AbstractMatrix{T}) where {T<:Integer}

Return a matrix whose columns generate the cone
intersection of the nonnegative orthant and the nullsapce of `N`.

# Examples
```jldoctest; setup = :(using CRNT)
julia> NonnegativeNullspaceCone([0 0; 0 0])
2×2 Array{Int64,2}:
 1  0
 0  1

julia> NonnegativeNullspaceCone([1 -1; 2 -2])
2×1 Array{Int64,2}:
 1
 1
```
"""
function NonnegativeNullspaceCone(N::AbstractMatrix{T}) where {T<:Integer}
    ## Migrating to Nemo
    Nnemo = Nemo.matrix(Nemo.FlintZZ,N)
    ## Computing nullsapce
    r, U = Nemo.nullspace_right_rational(Nnemo)
    ## Comming back to julia
    nullsp = T.(Array(U[:,1:r]))
    ## The vector space generated by nullsp as a cone
    rays1 = transpose(hcat(nullsp,-nullsp))
    c1 = Polymake.polytope.Cone(INPUT_RAYS=rays1)
    ## The nonnegative orthant as a cone
    rays2 = Matrix{T}(I, size(N,2), size(N,2))
    c2 = Polymake.polytope.Cone(INPUT_RAYS=rays2)
    ## Compute intersection of c1 c2 and save RAYS
    d = Polymake.polytope.intersection(c1,c2).RAYS
    ## Convert to integers
    d = Polymake.@convert_to Matrix{Integer} d
    ## Return the intersection as a Base.Array matrix of integers
    return T.(transpose(Array(d)))
end

include("TypesCompatibilites.jl")
include("NemoMatrices.jl")
include("StoichiometricMatrices.jl")

end
