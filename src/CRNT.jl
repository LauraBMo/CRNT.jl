module CRNT

## ENV["NEMO_PRINT_BANNER"] = false
using Nemo # matrix, FlintIntegerRing, nullspace
import AbstractAlgebra
import Polymake # polytope.cone,intersection
import LinearAlgebra: I, dot, Diagonal
import PolynomialRoots


# function __init__()
# end

###############################################################################
#                                  Extensions                                 #
###############################################################################
include("Base.jl")
include("Nemo.jl")
# include("Reduce.jl")

###############################################################################
#                                   Packages                                  #
###############################################################################

include("Cones.jl")
include("ConvexParameters.jl")
include("Stability.jl")
include("StoichiometricMatrices.jl")
include("CollectAndFinding.jl")
include("Maple.jl")

end
