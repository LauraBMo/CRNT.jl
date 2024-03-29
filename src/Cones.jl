
export cone_positiveorthant,
    cone_vectorspace,
    cone_positivenullspace,
    Newtonpolytope,
    verticesof,
    raysof,
    linearcombination,
    rays_outernormalcone,
    findnegativepoint,
    printfound,
    pRoots_qPossitive

## Nemo: matrix, FlintZZ, nullspace
## Polymake: polytope:cone,intersection
## LinearAlgebra: I

"""

    cone_positiveorthant(n)

Return the cone (Polymake big object) corresponding to the nonnegative orthant of R^n.

"""
function cone_positiveorthant(n)
    rays = cat(ones(Int, n)..., dims=(1, 2))
    return Polymake.polytope.Cone(INPUT_RAYS=rays)
end

"""

    cone_vectorspace(M::AbstractMatrix{T}) where {T <: Integer}

Return the cone (Polymake big object) corresponding to vector space generated by the columns of M.

"""
function cone_vectorspace(M::AbstractMatrix{T}) where {T <: Integer}
    rays = transpose(hcat(M, -M))
    return Polymake.polytope.Cone(INPUT_RAYS=rays)
end

"""

    cone_positivenullspace(N::AbstractMatrix{T}) where {T<:Integer}

Return the cone (Polymake big object) intersection of the nonnegative orthant and the nullspace of `N`.

"""
function cone_positivenullspace(N::AbstractMatrix{T}) where {T <: Integer}
    nullspace = Nemo.nullspace_right_rational(N)
    return Polymake.polytope.intersection(cone_positiveorthant(size(N, 2)), cone_vectorspace(nullspace))
end

function Newtonpolytope(p)
    points = transpose(collect_homoexponent_vectors(p))
    return Polymake.polytope.Polytope(POINTS=points)
end

function verticesof(polyt, polytname="")
    length(polytname) > 0 && print("Computing vertices of $(polytname), it may take some time.\n")
    return convert_to_array(polyt.VERTICES[:,2:end])
end

function raysof(cone, conename="")
    length(conename) > 0 && print("Computing rays of $(conename), it may take some time.\n")
    return Int.(integermultiple(Rational.(transpose(Array(cone.RAYS)))))
end

## Input:  polytope P, vertex v::Int.
## Output: point/vector in the outer normal cone
function linearcombination(A, coeffs=ones(Int, size(A, 2)))
    return A * coeffs
end

function rays_outernormalcone(polyt, vertex)
    cone = Polymake.polytope.normal_cone(polyt, vertex - 1, outer=1)
    return raysof(cone)
end

rays_outernormalcone(polyt) =  vertex -> rays_outernormalcone(polyt, vertex)

function findnegativepoint(p)
    Newtonp = Newtonpolytope(p)
    ver = verticesof(Newtonp, "Newton polytope")
    negver = Findallcols(isnegative, p, ver)
    points = []
    for v in negver
        rays = rays_outernormalcone(Newtonp, v)
        exponents = integermultiple(rays * ones(Int, size(rays, 2)))
        push!(points, collect_realpositiveroots(p, exponents))
    end
    return points
end

## Find roots of p for which q is possitive
function pRoots_qPossitive(p, q, nattemps::Integer=10, randbound::Integer=50)
    Newp = Newtonpolytope(p)
    Newq = Newtonpolytope(q)
    Vp = verticesof(Newp)
    Vq = verticesof(Newq, "q")
    println("Vertices of the Newton polytopes computed!!!!\n")
    #     P = NPq
    #     idx = posVq
    # else
    #     Cones = [Polymake.polytope.normal_cone(NPq, v-1, outer=1) for v in posVq]
    #     P = NPp
    #     idx = negVp
    #     ## Flip matrices of vertices to homogenize printing (the indices are flip)
    #     aux = VNPp
    #     VNPp = VNPq
    #     VNPq = aux
    #     end
    for i in posVq
        println("Computing cone to the $(i)th positive vertex of q\n")
        c1 = Polymake.polytope.normal_cone(NPq, i - 1, outer=1)
        for (nc2, c2) in enumerate(Conesp)
            rays = raysof(Polymake.polytope.intersection(c1, c2))
            r = size(rays, 1)
            if r > 0
                println("Intersecting cones found :)")
                found = false
                point = IntPointIncone(rays, ones(Int, 1, r))
                println("Computing real positive roots\n")
                proots = CollectallPossitiveRoots(p, point)
                qvals = [Float64(evaluate(q, (Nemo.fmpq).(pt))) for pt in proots]
                if findfirst(>(0), qvals) == nothing
                    j = 1;
                    while !found && j < nattemps
                        j += 1
                        println("Computing point in cone\n")
                        point = IntPointIncone(rays, rand(1:randbound, 1, r))
                        println("Computing real positive roots\n")
                        proots = CollectallPossitiveRoots(p, point)
                        qvals = [Float64(evaluate(q, (Nemo.fmpq).(pt))) for pt in proots]
                        found = findfirst(>(0), qvals) != nothing
                    end
                end
                if found
                    printfound(VNPp[negVp[nc2],:], VNPq[i,:], point, proots, qvals)
                end
            end
        end
    end
end


function printfound(Vp, Vq, point, proots, qvals)
    print("==============================================\n")
    print("==============================================\n")
    print("==============  Points found  ================\n")
    print("Negative Vertex p\n")
    print("$(Vp)\n")
    print("Possitive Vertex q\n")
    print("$(Vq)\n")
    print("Point in both outer normal cones\n")
    print("$(point)\n")
    for j in findall(>(0), qvals)
        print("$(proots[j])\n")
    end
end

function pRoots_qPossitive(p, NPp, VNPp, q, NPq, VNPq, nattemps::Integer=10, randbound::Integer=50)
    println("Finding negative/possitive vertices\n")
    isneg(x) = x < 0
    ispos(x) = 0 < x
    negVp = Findallrows(isneg, p, VNPp)
    posVq = Findallrows(ispos, q, VNPq)
    println("Computing cones to negative vertices of p\n")
    Conesp = [Polymake.polytope.normal_cone(NPp, v - 1, outer=1) for v in negVp]
    Listpt = []
    for i in posVq
        print("==============================================\n")
        print("==============================================\n")
        println("Computing the cone to the $(i)th positive vertex of q\n")
        c1 = Polymake.polytope.normal_cone(P, i - 1, outer=1)
        for (nc2, c2) in enumerate(Conesp)
            # print("==============================================\n")
            # print("==============================================\n")
            # println("Intersecting the $(nc2) first cone and the $i second\n")
            rays = raysof(Polymake.polytope.intersection(c1, c2))
            r = size(rays, 1)
            if r > 0
                # println("Intersecting cones found :)")
                found = false
                # println("Computing point in cone (coeffs 1)\n")
                point = IntPointIncone(rays, ones(Int, 1, r))
                # print("Point in both outer normal cones\n")
                # print("$(point)\n")
                # println("Computing real positive roots\n")
                proots = CollectallPossitiveRoots(p, point)
                qvals = [Float64(Nemo.evaluate(q, (Nemo.fmpq).(vec(pt)))) for pt in proots]
                pvals = [Float64(Nemo.evaluate(p, (Nemo.fmpq).(vec(pt)))) for pt in proots]
                found = findfirst(>(0), qvals) != nothing
                if !found
                    j = 1;
                    while !found && j < nattemps
                        j += 1
                        # println("Computing point in cone (rand coeffs, attempt $j)\n")
                        point = IntPointIncone(rays, rand(1:randbound, 1, r))
                        # print("Point in both outer normal cones\n")
                        # print("$(point)\n")
                        # println("Computing real positive roots\n")
                        proots = CollectallPossitiveRoots(p, point)
                        qvals = [Float64(Nemo.evaluate(q, (Nemo.fmpq).(vec(pt)))) for pt in proots]
                        pvals = [Float64(Nemo.evaluate(p, (Nemo.fmpq).(vec(pt)))) for pt in proots]
                        # print("Values of p for the roots\n")
                        # print("$(pvals)\n")
                        # print("Values of q for the roots\n")
                        # print("$(qvals)\n")
                        found = findfirst(>(0), qvals) != nothing
                    end
                end
                if found
                    printfound(VNPp[negVp[nc2],:], VNPq[i,:], point, proots, qvals)
                    Listpt = [Listpt; [[proots, qvals, pvals, point, VNPp[negVp[nc2],:], VNPq[i,:]]]]
                end
            end
        end
    end
    return Listpt
end
