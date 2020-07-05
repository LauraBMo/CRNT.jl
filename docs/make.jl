using CRNT
using Documenter

makedocs(;
    modules=[CRNT],
    authors="Laura Brustenga i Moncusí <brust@math.ku.dk> and contributors",
    repo="https://github.com/LauraBMo/CRNT.jl/blob/{commit}{path}#L{line}",
    sitename="CRNT.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://LauraBMo.github.io/CRNT.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/LauraBMo/CRNT.jl",
)
