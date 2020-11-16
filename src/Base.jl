
export convert_to_array,
    anynonzero,
    isrealof,
    isnegative,
    filter_complex,
    filter_negative,
    findfirstnonzero,
    findpivotsof,
    nonzeroslicesof,
    dropzeroslices,
    integermultiple


###############################################################################
#                             Types compatibilites                            #
###############################################################################

## Add float arithmetic compatibility with rational type
Base.Float64(x::Nemo.fmpq) = Base.Float64(Rational(x))
Base.Float32(x::Nemo.fmpq) = Base.Float32(Rational(x))
Base.Float16(x::Nemo.fmpq) = Base.Float16(Rational(x))
Base.BigFloat(x::Nemo.fmpq) = Base.BigFloat(Rational(x))

## In case we want to push! a fmpq coeff to a float array
Base.convert(::Type{T}, x::Nemo.fmpq) where {T <: AbstractFloat} = T(Rational(x))
Base.convert(::Type{T}, x::Nemo.fmpz) where {T <: AbstractFloat} = T(Int(x))

## From float to fmpq
(a::Nemo.FlintRationalField)(x::T) where {T <: AbstractFloat} = Nemo.fmpq(rationalize(x))
Nemo.fmpq(x::T) where {T <: AbstractFloat} = Nemo.fmpq(rationalize(x))

## From Nemo integers to Juila rationals
Rational{T}(x::Nemo.fmpz) where {T <: Integer} = Rational{T}(T(x))

## From Polymake matrices to Juila arrays
function convert_to_array(M::Polymake.Matrix, ::T=1) where {T}
    mat = Polymake.@convert_to Matrix{Integer} M
    return T.(transpose(Array(mat)))
end

## From Nemo matrices to Juila arrays
function convert_to_array(M::S, ::T=1) where {S <: MatElem,T}
    return T.(Array(M))
end


###############################################################################
#                         Miscelaneos small functions                         #
###############################################################################

# Nemo.isnonzero(x) = (x) != (zero(x))
# use !iszero instead (it is performance equivalent and more clean)
anynonzero(A) = any((!iszero).(A))

isrealof(rtol) = z -> abs(imag(z)) < rtol #

filter_complex(A, rtol) = real.(filter(isrealof(rtol), A))
filter_negative(A) = filter(isnegative, A)
# filter_positiveorthant(A, rtol) = filternegative(filtercomplex(A, rtol))

ispositive(x) = x>zero(x)
isnegative(x) = x<zero(x)

isnonnegative(x) = !(isnegative(x))
isnonpositive(x) = !(ispositive(x))

function getcoeff(term::Tuple{T,S}) where {T,S}
    return term[1]
end

function getexponent(term::Tuple{T,S}) where {T,S}
    return term[2]
end

function filterterms(predicate, p::MPolyElem)
	return Iterators.filter(predicate, dissect(p))
end

function findallcoeffs(predicate, p::MPolyElem)
    return findall(predicate, collect(coeffs(p)))
end

###############################################################################
#                             Array manipulations                             #
###############################################################################

function findfirstnonzero(A)
    return findfirst((!iszero), A)
end

function findpivotsof(W)
    # return findfirstnonzero.(eachrow(W)) ## it works just for matrices, code below is more generic.
    # return mapslices(findfirstnonzero, W, dims=collect(2:ndims(W)))
    return n -> findfirstnonzero.(eachslice(W, dims=n))
end

function nonzeroslicesof(A)
    # return n->vec(mapslices(anynonzero, M, dims=deleteat!(collect(1:ndims(M)),n)))
    return n -> anynonzero.(eachslice(A, dims=n))
end

function dropzeroslices(A)
    indeces = nonzeroslicesof(A).(1:ndims(A))
    return A[indeces...]
end

"""

    integermultiple(A)

Given an array of rational numbers 'A', returns a multiple 'λA' with integer entries, where 'λ' is the minimal integer with this property.
"""
function integermultiple(A)
    return Int.(abs(lcm(denominator.(A))) .* A)
end
