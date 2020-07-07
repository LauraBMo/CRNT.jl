
export StoichiometricCoeff, StoichiometricSources, StoichiometricTarget, Stoichiometric

function StoichiometricCoeff(net::AbstractMatrix, xs::AbstractVector; getcoefficient=DynamicPolynomials.coefficient)
    return hcat([getcoefficient(p,x) for x in xs, p in net[:,1]],
                [getcoefficient(p,x) for x in xs, p in net[:,2]])
end

function StoichiometricSources(Y::AbstractMatrix)
    return Y[:,1:div(end,2)]
end

function StoichiometricTarget(Y::AbstractMatrix)
    return Y[:,div(end,2)+1:end]
end

function Stoichiometric(Y::AbstractMatrix)
    return matYtarget(Y) - matYsources(Y)
end
