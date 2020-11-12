


ispositive(x) = x>zero(x)
isnegative(x) = x<zero(x)

isnonnegative(x) = !(isnegative(x))
isnonpositive(x) = !(ispositive(x))

Nemo.isnonzero(x) = !(iszero(x))

function findfirstnonzero(r)
    return findfirst(isnonzero,r)
end

function findpivots(W::AbstractMatrix)
    return findfirstnonzero.(eachrow(W))
end

function termsof(p::MPolyElem)
    return zip(coeffs(p),exponent_vectors(p))
end

function getcoeff(term::Tuple{T,S}) where {T,S}
    return term[1]
end

function getexponent(term::Tuple{T,S}) where {T,S}
    return term[2]
end

function filterterms(predicate, p::MPolyElem)
	return Iterators.filter(predicate, termsof(p))
end

function findallcoeffs(predicate, p::MPolyElem)
    return findall(predicate, collect(coeffs(p)))
end
