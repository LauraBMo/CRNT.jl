var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = CRNT","category":"page"},{"location":"#CRNT","page":"Home","title":"CRNT","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [CRNT]","category":"page"},{"location":"#AbstractAlgebra.Generic.matrix-Tuple{AbstractAlgebra.Ring,AbstractArray{T,1} where T}","page":"Home","title":"AbstractAlgebra.Generic.matrix","text":"matrix(R, V::AbstractVector{T}) where {T}\n\nGiven a Julia vector V of entries, construct the corresponding AbstractAlgebra.jl one-column matrix over the given ring R, assuming all the entries can be coerced into R.\n\n\n\n\n\n","category":"method"},{"location":"#CRNT.Diagonal-Union{Tuple{T}, Tuple{AbstractAlgebra.Ring,AbstractArray{T,1}}} where T","page":"Home","title":"CRNT.Diagonal","text":"Diagonal(R, V::AbstractVector{T}) where {T}\n\nGiven a Julia vector V of entries, construct the corresponding AbstractAlgebra.jl diagonal matrix over the given ring R, assuming all the entries can be coerced into R.\n\n\n\n\n\n","category":"method"},{"location":"#CRNT.Jacobian-Union{Tuple{T}, Tuple{AbstractAlgebra.Ring,AbstractAlgebra.MatElem{T}}, Tuple{AbstractAlgebra.Ring,AbstractAlgebra.MatElem{T},Any}} where T<:AbstractAlgebra.RingElem","page":"Home","title":"CRNT.Jacobian","text":"Jacobian(R, M)\nJacobian(R, M, vars=1:length(Nemo.gens(R)))\n\nReturns the Jacobian matrix of a one-column matrix of polynomials M with respect to the generators of R indexed by vars. When vars is omitted all the generators of R are used.\n\nExamples\n\njulia> using Nemo\n\njulia> R, vars = PolynomialRing(ZZ, vcat([\"k$i\" for i in 1:5], [\"x$i\" for i in 1:4]));\n\njulia> M = matrix(R, [vars[1]*vars[2]*vars[6]-vars[8]; vars[3]*vars[9]+2*vars[7]])\n[k1*k2*x1-x3]\n[ k3*x4+2*x2]\n\njulia> Jacobian(R, M, 6:9)\n[k1*k2  0  -1   0]\n[    0  2   0  k3]\n\njulia> Jacobian(R,M)\n[k2*x1  k1*x1   0  0  0  k1*k2  0  -1   0]\n[    0      0  x4  0  0      0  2   0  k3]\n\n\n\n\n\n","category":"method"},{"location":"#CRNT.NonnegativeNullspaceCone-Union{Tuple{AbstractArray{T,2}}, Tuple{T}} where T<:Integer","page":"Home","title":"CRNT.NonnegativeNullspaceCone","text":"NonnegativeNullspaceCone(N::AbstractMatrix{T}) where {T<:Integer}\n\nReturn a matrix whose columns generate the cone intersection of the nonnegative orthant and the nullsapce of N.\n\nExamples\n\njulia> NonnegativeNullspaceCone([0 0; 0 0])\n2×2 Array{Int64,2}:\n 1  0\n 0  1\n\njulia> NonnegativeNullspaceCone([1 -1; 2 -2])\n2×1 Array{Int64,2}:\n 1\n 1\n\n\n\n\n\n","category":"method"}]
}
