
export
    StoichiometricCoeffs,
    StoichiometricSources,
    StoichiometricTargets,
    StoichiometricMatrix,
    Velocity,
    Velocities,
    VectorVelocities,
    Speciesformationrate,
    SpeciesformationTotalRing,
    SpeciesformationKsRing,
    SpeciesformationRings,
    trivialConservativeLaws

# Matrix Y
function StoichiometricCoeffs(net::AbstractMatrix, xs::AbstractVector, getcoefficient=Nemo.coeff)
    return Int.(hcat([getcoefficient(p,x) for x in xs, p in net[:,1]],
                     [getcoefficient(p,x) for x in xs, p in net[:,2]]))
end

function StoichiometricSources(Y::AbstractMatrix)
    return Y[:,1:div(end,2)]
end

function StoichiometricTargets(Y::AbstractMatrix)
    return Y[:,div(end,2)+1:end]
end

# Matrix N
function StoichiometricMatrix(Y::AbstractMatrix)
    return StoichiometricTargets(Y) - StoichiometricSources(Y)
end

##########################
##########################
##########################
######## W
########

function trivialConservativeLaws(N::AbstractMatrix{T}) where {T<:Integer}
    nTs, W = left_kernel(matrix(FlintIntegerRing(),N))
    return nTs, Int.(Array(hnf(W)))
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
function Velocity(R::MPolyRing, v::AbstractVector)
    newV = MPolyBuildCtx(R)
    push_term!(newV, one(R), v)
    return finish(newV)
end

function Velocities(R::MPolyRing)
    return v->Velocity(R, v)
end

function columns(A::AbstractMatrix)
    return [A[:, i] for i in 1:size(A,2)]
end

function VectorVelocities(R::MPolyRing, Y::AbstractMatrix)
    return Velocities(R).(columns(Y))
end

function SpeciesformationTotalRing(net, nxs::Integer)
    nks = size(net, 1)
    return PolynomialRing(FlintRationalField(), vcat(["x$i" for i in 1:nxs], ["k$i" for i in 1:nks]))
end

function SpeciesformationKsRing(net)
    nks = size(net, 1)
    return PolynomialRing(FlintRationalField(), ["k$i" for i in 1:nks])
end

function SpeciesformationRings(net, nxs::Integer)
    return (SpeciesformationKsRing(net)..., PolynomialRing(SpeciesformationKsRing(net)[1], ["x$i" for i in 1:nxs])...)
end

function Speciesformationrate(Rx::MPolyRing, Rk::MPolyRing,
                              N::AbstractMatrix, Y::AbstractMatrix)
    V = VectorVelocities(Rx, Y)
    K = Diagonal(gens(Rk))
    return N*K*V
end

isnonzero(x) = (x)!=(zero(typeof(x)))

function findfirstnonzero(r)
    return findfirst(isnonzero,r)
end

function findpivots(W::AbstractMatrix)
    return findfirstnonzero.(eachrow(W))
end

function SpeciesformationrateInStclass!(f::AbstractVector,
                                        R::MPolyRing,
                                        W::AbstractMatrix,
                                        xs::UnitRange{T}=1:size(W,2),
                                        ts::UnitRange{T}=.+(length(xs),1:size(W,1))) where {T<:Integer}
    Wx = W*gens(R)[xs] - gens(R)[ts]
    for (i,p) in enumerate(findpivots(W))
        f[p] = Wx[i]
    end
    return f
end

function SpeciesformationrateInStclass(f::AbstractVector,
                                       R::MPolyRing,
                                       W::AbstractMatrix,
                                       xs::UnitRange{T}=1:size(W,2),
                                       ts::UnitRange{T}=.+(length(xs),1:size(W,1))) where {T<:Integer}
    F = f
    return SpeciesformationrateInStclass!(F, R, W, xs, ts)
end

function SystemF(net, nxs)
    Y = StoichiometricCoeffs(net, nxs)
    Ys = StoichiometricSources(Y)
    N = StoichiometricMatrix(Y)
