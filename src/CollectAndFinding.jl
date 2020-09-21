
export
    homogenize,
    homoexponent_vectors,
    CollectallPossitiveRoots,
    CollectHomoiterInmatrows,
    Collectallrows

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
    return pushfirst!(Vector(n[1]), oneunit(Base.eltype(n[1]))), n[2]
end

Base.@propagate_inbounds function Base.iterate(h::Homogenize)
    n = iterate(h.itr)
    n === nothing && return n
    return pushfirst!(Vector(n[1]), oneunit(Base.eltype(n[1]))), n[2]
end



"""
    homoexponent_vectors(p::MPolyElem)

An iterator for the homogenized (with an extra first component set to '1' exponent vectors of the multivariate polynomial `p`.
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

julia> collect(homoexponent_vectors(p))
2-element Array{Any,1}:
 [1, 1, 1, 0, 0, 0, 1, 0, 0, 0]
 [1, 0, 0, 0, 0, 0, 0, 0, 1, 0]
```
"""
homoexponent_vectors(p::MPolyElem) = homogenize(exponent_vectors(p))

function collect_homoiter(iter)
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
        v[2:end,state-1] = u
        next = iterate(iter, state)
    end
    return transpose(v)
end

"""
    CollectExponent_homovectors(p::MPolyElem)

Collects the homogenized (with an extra first component set to 1) exponents of `p` in the rows of a matrix.
The returned matrix can be used to compute the Newton polytope of `p`.
# Examples
```jldoctest; setup = :(using CRNT, Nemo)
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> q = vars[1]*vars[2]*vars[6]-vars[8]
k1*k2*x1-x3

julia> collect(exponent_vectors(q))
2-element Array{Array{Int64,1},1}:
 [1, 1, 0, 0, 0, 1, 0, 0, 0]
 [0, 0, 0, 0, 0, 0, 0, 1, 0]

julia> CollectExponent_homovectors(q)
2×10 Transpose{Int64,Array{Int64,2}}:
 1  1  1  0  0  0  1  0  0  0
 1  0  0  0  0  0  0  0  1  0
```

Recall it is approximately ten times faster than the equivalent
    transpose(reduce(hcat,collect(exponent_homovectors(q))))
"""
function CollectExponent_homovectors(p)
    ## Save the iterator (this function just works for iters whose states are Int and start at ZERO)
    iter = exponent_vectors(p)
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

"""
    Collectallrows(predicate, p::MPolyElem, V::AbstractMatrix)

Returns an `Array` of all the exponents of terms of `p` for which the exponent is a row of `V` and the coefficient satisfies `predicate`.

# Examples
```jldoctest; setup = :(using CRNT)
julia> using Nemo

julia> R, vars = PolynomialRing(ZZ, vcat(["k\$i" for i in 1:5], ["x\$i" for i in 1:4]));

julia> q = vars[1]*vars[2]*vars[6]-vars[3]*vars[8]+vars[5]*vars[9]-2*vars[1]*vars[3]*vars[5]^2
k1*k2*x1-2*k1*k3*k5^2-k3*x3+k5*x4

julia> collect(exponent_vectors(q))
4-element Array{Array{Int64,1},1}:
 [1, 1, 0, 0, 0, 1, 0, 0, 0]
 [1, 0, 1, 0, 2, 0, 0, 0, 0]
 [0, 0, 1, 0, 0, 0, 0, 1, 0]
 [0, 0, 0, 0, 1, 0, 0, 0, 1]

julia> V =  [[1 1 0 0 0 1 0 0 0]; [0 0 1 0 0 0 0 1 0]; [1 0 1 0 2 0 0 0 0]];

julia> Collectallrows(x->Nemo.isless(x,0), q, V)
2-element Array{Array{Int64,1},1}:
 [1, 0, 1, 0, 2, 0, 0, 0, 0]
 [0, 0, 1, 0, 0, 0, 0, 1, 0]
```
"""
function Collectallrows(predicate, p::MPolyElem, V::AbstractMatrix)
    predisinV(x) = predicate(x[1]) && x[2] in eachrow(V)
    getexponent(x) = x[2]
    return getexponent.(Iterators.filter(predisinV, zip(coeffs(p),exponent_vectors(p))))
end

function Findallrows(predicate, p::MPolyElem, V::AbstractMatrix)
    perdexpposs = findall(predicate, collect(coeffs(p)))
    perdexp = [iterate(exponent_vectors(p),n-1)[1] for n in perdexpposs]
    return findall(in(perdexp), collect(eachrow(V)))
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

function CollectallRootspositiveorthant(p, inNormalCone, realtol::Real=1e-7)
    ## Compute exponent of t (to know the minimum in case it is negative).
    texponents = [dot(inNormalCone,exp) for exp in exponent_vectors(p)]
    maxdeg = maximum(texponents)
    mindeg = minimum(texponents)
    ## If mindeg>0 we divide by t^(mindeg), otherwise we multiply.
    ## Both correspond to shift by -mindeg!
    poly = zeros(maxdeg-mindeg+1)
    for (c,i) in zip(coeffs(p),texponents)
        poly[i-mindeg+1] = c
    end
    # mindeg = findfirst(!=(0), poly)
    # roots = PolynomialRoots.roots(BigFloat.(poly[mindeg:end]))
    # roots = PolynomialRoots.roots(poly[mindeg:end])
    roots = PolynomialRoots.roots(BigFloat.(poly))
    realposs = filter(>(0), real.(filter(x->(abs(imag(x))<realtol), roots)))
    return [t.^inNormalCone for t in realposs]
end
