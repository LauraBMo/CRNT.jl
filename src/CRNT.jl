module CRNT

## ENV["NEMO_PRINT_BANNER"] = false
using Nemo # matrix, FlintIntegerRing, nullspace
import AbstractAlgebra: Ring
import Polymake # polytope.cone,intersection
import LinearAlgebra: I, dot
import PolynomialRoots

# function __init__()
# end


include("Cones.jl")
include("ConvexParameters.jl")
include("Stability.jl")
include("TypesCompatibilites.jl")
include("NemoMatrices.jl")
include("StoichiometricMatrices.jl")
include("CollectAndFinding.jl")

end
