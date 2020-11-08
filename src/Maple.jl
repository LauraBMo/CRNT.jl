

function toMaple(net, nxs, file::String)
    S = stoichiometriccoeffs(net, nxs) ## S may have zero rows and cols
    N = stoichiometricmatrix(S) ## Has no zero row or col
    nts, W = conservativelaws(N)
    Y = kineticorder(S)
    E = cone_positivenullspace(N)
    open(file, "w") do io
        write(io, "read(\"ModelsMatricies/Procs.mpl\"):\n")
        @matrixtoMaple io Y
        @matrixtoMaple io N
        @matrixtoMaple io W
        @matrixtoMaple io E
        xstoMaple(io, stoichiometricsources(S))
        kstoMaple(io, stoichiometricsources(S))
        systemtoMaple(io)
        WsystemtoMaple(io, nts, W)
    end
end

function matrixtoMaple(io, M, name)
    write(io, "\n\n$(name) := Matrix$(size(M)):\n\n")
    for i in CartesianIndices(M)
        if M[i] != 0
            write(io, "$(name)[$(i[1]),$(i[2])] := $(M[i]):\n")
        end
    end
end

macro matrixtoMaple(io, M)
    str = string(M)
    return quote
        matrixtoMaple($io, $M, $str)
    end
end

function xstoMaple(io, stoichiometricsources)
    xs = (1:(size(stoichiometricsources, 1)))[nonzerorows(stoichiometricsources)]
    write(io, "\n\nnxs := $(size(xs, 1)):\n")
    write(io, "xs := [seq(x[i], i = [")
    for x in xs[1:(end - 1)]
        write(io, "$x, ")
    end
    write(io, "$(xs[end])])];\n\n")
    write(io, "depenvars := [seq(x[i], i = [")
    for x in xs[2:(end - 2)]
        write(io, "$x, ")
    end
    write(io, "$(xs[end - 1])])];\n\n")
end

function kstoMaple(io, stoichiometricsources)
    ks = (1:(size(stoichiometricsources, 2)))[nonzerocols(stoichiometricsources)]
    write(io, "\n\nnks := $(size(ks, 1)):\n")
    write(io, "ks := [seq(k[i], i = [")
    for k in ks[1:(end - 1)]
        write(io, "$k, ")
    end
    write(io, "$(ks[end])])];\n\n")
end

function systemtoMaple(io)
    write(io, "v := Velocities(Y,xs):\n")
    write(io, "digK := Matrix(nks):\n")
    write(io, "for i to nks do digK[i,i] := ks[i] end do:\n")
    write(io, "S := N.digK.v:\n")
    write(io, "Seq:= equfy(S);\n")
end

function WsystemtoMaple(io, nts, W)
    write(io, "\n\nSw := copy(S):\n")
    write(io, "Wx := (W.(Vector[column](xs))) - Vector[column]([seq(T[i], i = 1 .. $(nts))]):\n")
    for (i, p) in enumerate(findpivots(W))
        write(io, "Sw[$p] := Wx[$i]:\n")
    end
    write(io, "Sweq:= equfy(Sw):\n")
    write(io, "J := VectorCalculus[Jacobian](Sw, xs):\n")
    write(io, "DJ := (-1)^(Rank(N))*Determinant(J):\n")
end
