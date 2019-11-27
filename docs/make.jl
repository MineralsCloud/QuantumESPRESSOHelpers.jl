using Documenter, QuantumESPRESSOHelpers

makedocs(;
    modules=[QuantumESPRESSOHelpers],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/MineralsCloud/QuantumESPRESSOHelpers.jl/blob/{commit}{path}#L{line}",
    sitename="QuantumESPRESSOHelpers.jl",
    authors="Qi Zhang <singularitti@outlook.com>",
    assets=String[],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOHelpers.jl",
)
