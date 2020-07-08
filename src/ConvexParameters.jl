
export
    JacobianConvexparameters,
    CharpolyCoeff,
    JacobianDeterminantConvexparameters

function JacobianConvexparameters(S::Ring,
                                  N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix,
                                  lpos::UnitRange{T} = 1:size(E,2),
                                  hpos::UnitRange{T} = .+(length(lpos),1:size(N,1))) where {T<:Integer}
    ## Migrating to AbstractAlgebra
    Nnemo = Nemo.matrix(S, N)
    Ynemo = Nemo.matrix(S, Y)
    ## Square matrix with h's diagonal
    H = Diagonal(S, Nemo.gens(S)[hpos])
    ## Diagonal with coordinates of E
    D = Diagonal(S, E*Nemo.gens(S)[lpos])
    ## Return the jacobian under convex parameters
    return Nnemo*D*Ynemo*H
end

function CharpolyCoeff(R::Ring, M::MatElem, ncoeff::Integer)
    P, x = PolynomialRing(R, "x")
    # q = coeff(divexact(charpoly(P,Jhl),x^mindeg),0);
    return coeff(charpoly(P, M), ncoeff)
end

function JacobianDeterminantConvexparameters(R::Ring,
                                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix,
                                             lpos::UnitRange{T} = 1:size(E,2),
                                             hpos::UnitRange{T} = .+(length(lpos),1:size(N,1)),
                                             mindeg::Integer = size(N,1)-rank(Nemo.matrix(R, N))) where {T<:Integer}
    return CharpolyCoeff(R, JacobianConvexparameters(R, N, Y, E, lpos, hpos), mindeg)
end
