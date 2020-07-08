
export
    StoichiometricCoeff,
    StoichiometricSources,
    StoichiometricTargets,
    StoichiometricMat,
    VectorVelocities,
    Speciesformationrate

function StoichiometricCoeff(net::AbstractMatrix, xs::AbstractVector; getcoefficient=DynamicPolynomials.coefficient)
    return hcat([getcoefficient(p,x) for x in xs, p in net[:,1]],
                [getcoefficient(p,x) for x in xs, p in net[:,2]])
end

function StoichiometricSources(Y::AbstractMatrix)
    return Y[:,1:div(end,2)]
end

function StoichiometricTargets(Y::AbstractMatrix)
    return Y[:,div(end,2)+1:end]
end

function StoichiometricMat(Y::AbstractMatrix)
    return StoichiometricTargets(Y) - StoichiometricSources(Y)
end

##############################################
##############################################
##############################################
############ Building matrices in Nemo for CRN

################### With k and x
## Velocities
## Input: Array with 'x' variables xs, and the kinetics matrix
## Output: Vector of velocities of the reaction network
## VectorVelocities(gens(R), Y)
function VectorVelocities(R::Ring, Y::AbstractMatrix,
                          xs::UnitRange{T}=1:size(Y,1)) where {T<:Integer}
    V = [prod(Nemo.gens(R)[xs].^c) for c in eachcol(Y)]
    return Nemo.matrix(R, V)
end

function Speciesformationrate(R::Ring, K::Ring,
                              N::AbstractMatrix, Y::AbstractMatrix,
                              xs::UnitRange{T}=1:size(Y,1), ks::UnitRange{T}=1:size(Y,2)) where {T<:Integer}
    # r = size(Y, 2)
    ## N in Nemo
    Nnemo = Nemo.matrix(R, N)
    ## Velocities
    V = VectorVelocities(R, Y, xs)
    ## Rates
    K = Diagonal(R, Nemo.gens(K)[ks])
    ## Finally the system
    return Nnemo*K*V
end
