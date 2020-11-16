
export JacobianConvexparameters,
    CharpolyCoeff,
    JacobianDeterminantConvexparameters

function JacobianConvexparameters(S::AbstractAlgebra.Ring,
                                  N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix,
                                  lpos::UnitRange{T}=1:size(E, 2),
                                  hpos::UnitRange{T}=.+(length(lpos), 1:size(N, 1))) where {T <: Integer}
    ## Migrating to AbstractAlgebra
    Nnemo = Nemo.matrix(S, N)
    Ynemo = Nemo.matrix(S, transpose(Y))
    ## Square matrix with h's diagonal
    H = Diagonal(S, Nemo.gens(S)[hpos])
    ## Diagonal with coordinates of E
    D = Diagonal(S, E * Nemo.gens(S)[lpos])
    ## Return the jacobian under convex parameters
    return Nnemo * D * Ynemo * H
end

function CharpolyCoeff(R::AbstractAlgebra.Ring, M::MatElem, ncoeff::Integer)
    P, x = Nemo.PolynomialRing(R, "x")
    # q = coeff(divexact(charpoly(P,Jhl),x^mindeg),0);
    return Nemo.coeff(charpoly(P, M), ncoeff)
end

function JacobianDeterminantConvexparameters(R::AbstractAlgebra.Ring,
                                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix,
                                             lpos::UnitRange{T}=1:size(E, 2),
                                             hpos::UnitRange{T}=.+(length(lpos), 1:size(N, 1)),
                                             mindeg::Integer=size(N, 1) - rank(Nemo.matrix(R, N))) where {T <: Integer}
    return CharpolyCoeff(R, JacobianConvexparameters(R, N, Y, E, lpos, hpos), mindeg)
end

function JacobianConvexparameters(H::AbstractAlgebra.Ring,
                                  L::AbstractAlgebra.Ring,
                                  N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix)
    ## Migrating to AbstractAlgebra
    Nnemo = Nemo.matrix(H, N)
    Ynemo = Nemo.matrix(H, transpose(Y))
    ## Square matrix with h's diagonal
    digH = Diagonal(H, Nemo.gens(H))
    ## Diagonal with coordinates of E
    digL = Diagonal(H, E * Nemo.gens(L))
    ## Return the jacobian under convex parameters
    return Nnemo * digL * Ynemo * digH
end

function JacobianDeterminantConvexparameters(H::AbstractAlgebra.Ring,
                                             L::AbstractAlgebra.Ring,
                                             N::AbstractMatrix, Y::AbstractMatrix, E::AbstractMatrix,
                                             mindeg::Integer=size(N, 1) - rank(Nemo.matrix(H, N)))
    return CharpolyCoeff(H, JacobianConvexparameters(H, L, N, Y, E), mindeg)
end
