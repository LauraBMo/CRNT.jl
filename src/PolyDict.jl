


function extractexponentsofterm(T::RExpr)
    return v->Algebra.deg(T,v)
end

function listofterms(P::RExpr)
    return RExpr.(split(String(P),r"\+|\-"))
end

function exponentsofpoly(P::RExpr, vars)
    return [extractexponentsofterm(t).(vars) for t in listofterms(P)]
end

function VectorOfCoeff(P::RExpr, var)
    return collect(RExpr.(parse(Algebra.coeff(P,var))))
end

# function VectorOfCoeff(P::RExpr, var)
#     return RExpr.(split(String(Algebra.coeff(q,var))[2:end-1], ","))
# end

function degOfPoly(P::RExpr)
    return v-> parse(Algebra.deg(P,v))
end

function PolyToTensor(P::RExpr,vars)
    dims = Tuple(1 .+ degOfPoly(P).(vars))
    N = size(dims,1)-1
    out = zeros(RExpr,dims)
    out[:,fill(1,N)...] = VectorOfCoeff(P,vars[1])
    for n in 1:N
        for I in CartesianIndices(dims[1:n])
            l = VectorOfCoeff(out[(Tuple(I)..., fill(1, N-n+1)...)...],vars[n+1])
            out[(Tuple(I)..., 1:size(l,1), fill(1, N-n)...)...] = l
        end
    end
    return out
end

function coeffOfcoeff(L::AbstractArray, var, d)
    Lnew = zeros(RExpr,(size(L)...,d))
    for i in CartesianIndices(L)
        l = listofcoeff(L[i],var)
        Lnew[(Tuple(i)...,1:size(l,1))...] = l
    end
    return Lnew
end

struct Rpoly
        P::RExpr
        L
        vars
end

function Base.iterate(P::Rpoly)
    return Rpoly(P.P,listofcoeff(P.P,P.vars[1]),P.vars[2:end]), 2
end

function coeffOfcoeff(L::AbstractArray, var, d)
    Lnew = zeros(RExpr,(size(L)...,d))
    for i in CartesianIndices(L)
        l = listofcoeff(L[i],var)
        Lnew[(Tuple(i)...,1:size(l,1))...] = l
    end
    return Lnew
end

function Base.iterate(P::Rpoly, state)
        if size(P.vars,1) > 0
            ## P.L is an Array of dim state-1 (state=2,vector;3,matrix;...)
            dnew = parse(Algebra.deg(Algebra.lterm(P.P,P.vars[1]),P.vars[1]))+1
            Lnew = coeffOfcoeff(P.L, P.vars[1], dnew)
            if size(P.vars,1) == 1
                Lnew = parse.(Lnew)
            end
            return Rpoly(P.P,Lnew,P.vars[2:end]), state+1
        else
        return nothing, state+1
    end
end

function PolyToTensor(P::RExpr,vars)
    next, state = iterate(Rpoly(P,[],vars))
    previous = next
    while next !== nothing
        previous = next
        next, state = iterate(next,state)
    end
    return previous.L
end
