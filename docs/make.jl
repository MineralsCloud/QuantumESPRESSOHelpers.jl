using QuantumESPRESSOHelpers
using Documenter

makedocs(;
    modules=[QuantumESPRESSOHelpers],
    authors="Qi Zhang <singularitti@outlook.com>",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOHelpers.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSOHelpers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOHelpers.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOHelpers.jl",
)
