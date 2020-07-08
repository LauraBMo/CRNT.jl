
export
    homogenize,
    exponent_homovectors

## Code adapted from "enumarate(iter)"
struct Homogenize{I}
    itr::I
end

"""
    homogenize(iter)

Given an iterator `iter` whose values are 1d arrays,
returns an iterator that yields `[1; x]` where `1` is a unit of the
element type of `x`, and `x` is a value of `iter`.

# Examples
```jldoctest; setup = :(using CRNT)
julia> x = [[1, 2], [3, 4], [5, 6]];

julia> collect(homogenize(x))
3-element Array{Array{Int64,1},1}:
 [1, 1, 2]
 [1, 3, 4]
 [1, 5, 6]
```
"""
homogenize(iter) = Homogenize(iter)

Base.length(h::Homogenize) = Base.length(h.itr)
Base.size(h::Homogenize) = Base.size(h.itr)
Base.eltype(::Type{Homogenize{I}}) where {I} = eltype(I)
Base.IteratorSize(::Type{Homogenize{I}}) where {I} = Base.IteratorSize(I)
Base.IteratorEltype(::Type{Homogenize{I}}) where {I} = Base.IteratorEltype(I)

Base.@propagate_inbounds function Base.iterate(h::Homogenize, state)
    n = iterate(h.itr, state...)
    n === nothing && return n
    return [oneunit(Base.eltype(n[1])); @view(n[1][:])], n[2]
end

Base.@propagate_inbounds function Base.iterate(h::Homogenize)
    n = iterate(h.itr)
    n === nothing && return n
    return [oneunit(Base.eltype(n[1])); @view(n[1][:])], n[2]
end



"""
    exponent_homovectors(p::MPolyElem)

An iterator for the homogenized (with an extra first component set to 1) exponent vectors of the multivariate polynomial `p`.
To retrieve an array, use `collect(exponent_homovectors(p))`.
# Examples
```jldoctest; setup = :(using CRNT, Nemo)
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> p = vars[1]*vars[2]*vars[6]-vars[8]
k1*k2*x1-x3

julia> collect(exponent_vectors(p))
2-element Array{Array{Int64,1},1}:
 [1, 1, 0, 0, 0, 1, 0, 0, 0]
 [0, 0, 0, 0, 0, 0, 0, 1, 0]

julia> collect(exponent_homovectors(p))
2-element Array{Any,1}:
 [1, 1, 1, 0, 0, 0, 1, 0, 0, 0]
 [1, 0, 0, 0, 0, 0, 0, 0, 1, 0]
```
"""
exponent_homovectors(p::MPolyElem) = homogenize(exponent_vectors(p))

"""
    CollectHomoiterInmatrows(iter)

Given an iterator `iter` whose values are 1d arrays,
collect its homogenized (with an extra first component set to 1) values in the rows of a matrix.
# Examples
```jldoctest; setup = :(using CRNT, Nemo)
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> p = vars[1]*vars[2]*vars[6]-vars[8]
k1*k2*x1-x3

julia> collect(exponent_vectors(p))
2-element Array{Array{Int64,1},1}:
 [1, 1, 0, 0, 0, 1, 0, 0, 0]
 [0, 0, 0, 0, 0, 0, 0, 1, 0]

julia> CollectHomoiterInmatrows(exponent_vectors(p))
2×10 Transpose{Int64,Array{Int64,2}}:
 1  1  1  0  0  0  1  0  0  0
 1  0  0  0  0  0  0  0  1  0
```

Recall it is approximately ten times more efficient than the equivalent
    transpose(reduce(hcat,collect(exponent_homovectors(q))))
"""
function CollectHomoiterInmatrows(iter)
    ## Save the iterator
    # iter = exponent_vectors(p)
    ## First iteration outside the loop to prealocate v
    (u, state) = iterate(iter)
    ## Prealocate v
    v = ones(eltype(u), size(u,1)+1, length(iter))
    # v = Matrix{eltype(u)}(undef, size(u,1)+1, length(iter))
    # v[1,:] .= one(eltype(u))
    v[2:end,1] = u
    next = iterate(iter, state)
    ## Start the loop
    while next !== nothing
        (u, state) = next
        v[2:end,state] = u
        next = iterate(iter, state)
    end
    return transpose(v)
end

getexponent(x) = x[2]
# predicate=(x->isnegative(x[1]))
function CollectallVertices(predicate, p, V)
    isnegvertex(x) = predicate(x) && x[2] in eachrow(V)
    return getexponent.(Iterators.filter(isnegvertex, zip(coeffs(p),exponent_vectors(p))))
end

# predicate=isnegative
function FindallVertices(predicate, p::MPolyElem, V::AbstractMatrix)
    negexpposs = findall(predicate, collect(coeffs(p)))
    negexp = [iterate(exponent_vectors(p),n-1)[1] for n in negexponentposs]
    return findall(in(negexp), collect(eachrow(V)))
end

# getcoeff(x) = x[1]
# isnegative(x) = x<0
# iscoeffnegative(x) = isnegative(getcoeff(x))
# negcoeffs_exponent(p) = Iterators.Filter(iscoeffnegative, zip(coeffs(p),exponent_vectors(p)))
# negexp = getexponent.(negcoeffs_exponent(p))

## Method to obtain negexp 4 times faster for the non-commented function
## (you need to define the polynomial q)
# using BenchmarkTools
# @benchmark getexponent.(negcoeffs_exponent(q))
# @benchmark [iterate(exponent_vectors(q),n-1)[1] for n in findall(isnegative, collect(coeffs(q)))]
