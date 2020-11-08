
@inline tuplejoin(x) = x
@inline tuplejoin(x, y) = (x..., y...)
@inline tuplejoin(x, y, z...) = (x..., tuplejoin(y, z...)...)



function exponentofRterm(T::RExpr)
    return v->parse(Algebra.deg(T,v))
end

function listofRterms(P::RExpr)
    return RExpr.(split(String(P),r"\+|\-"))
end

function exponentofRpoly(P::RExpr, vars)
    return [exponentofRterm(t).(vars) for t in listofRterms(P)]
end

function Rcoeffs(P::RExpr, var)
    return collect(RExpr.(parse(Algebra.coeff(P,var))))
end

# function Rcoeffs(P::RExpr, var)
#     return RExpr.(split(String(Algebra.coeff(q,var))[2:end-1], ","))
# end

function degOfPoly(P::RExpr)
    return v -> parse(Algebra.deg(P,v))
end

function RpolyToTensor(P::RExpr,vars)
    dims = Tuple(1 .+ degOfPoly(P).(vars))
    N = size(dims,1)-1
    out = zeros(RExpr,dims)
    out[:,fill(1,N)...] = Rcoeffs(P,vars[1])
    for n in 1:N
        for I in CartesianIndices(dims[1:n])
            # l = Rcoeffs(out[(Tuple(I)..., fill(1, N-n+1)...)...],vars[n+1])
            l = Rcoeffs(out[I, fill(1, N-n+1)...],vars[n+1])
            out[I, 1:size(l,1), fill(1, N-n)...] = l
        end
    end
    return out
end

function RpolyToNemo(RP::RExpr,Rvars,NR::MPolyRing,Nvars=Nemo.gens(R))
    exp = exponentofRpoly(RP, Rvars)
    coeffs = parse.(RpolyToTensor(RP, Rvars))
    newP = AbstractAlgebra.MPolyBuildCtx(NR)
    for e in exp
        AbstractAlgebra.push_term!(newP, base_ring(NR)(coeffs[(1 .+ e)...]), e)
    end
    return AbstractAlgebra.finish(newP)
end

function NemoToRpoly(NP)
    return RExpr(NP)
end
