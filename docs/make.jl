using QuantumESPRESSOHelpers
using Documenter

DocMeta.setdocmeta!(QuantumESPRESSOHelpers, :DocTestSetup, :(using QuantumESPRESSOHelpers); recursive=true)

makedocs(;
    modules=[QuantumESPRESSOHelpers],
    authors="singularitti <singularitti@outlook.com> and contributors",
    repo="https://github.com/MineralsCloud/QuantumESPRESSOHelpers.jl/blob/{commit}{path}#{line}",
    sitename="QuantumESPRESSOHelpers.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://MineralsCloud.github.io/QuantumESPRESSOHelpers.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MineralsCloud/QuantumESPRESSOHelpers.jl",
)
