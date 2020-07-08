module CRNT

using Nemo # matrix, FlintIntegerRing, nullspace
import AbstractAlgebra: Ring
import Polymake # polytope.cone,intersection
import LinearAlgebra: I, dot
import PolynomialRoots: roots

include("Cones.jl")
include("TypesCompatibilites.jl")
include("NemoMatrices.jl")
include("StoichiometricMatrices.jl")

end
