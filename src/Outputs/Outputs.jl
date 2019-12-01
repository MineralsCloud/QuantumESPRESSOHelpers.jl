module Outputs

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: GREEN_FG
using QuantumESPRESSOBase.Cards.PWscf: AtomicPositionsCard
using QuantumESPRESSOParsers.Outputs.PWscf

export PWOutput, output_parser

struct PWOutput end

function output_parser(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:PWOutput}
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    print(terminal, GREEN_FG("Please give the absolute path to your output file: "))
    path = abspath(strip(readline(terminal)))
    str = try
        open(path, "r") do io
            read(io, String)
        end
    catch
        @warn("File '$path' not found!")
    end
    calculation = calculations[request(
        terminal,
        GREEN_FG("What exact calculation is this output?") |> string,
        RadioMenu(collect(values(calculations))),
    )]
    if calculation ∈ ("relax", "vc-relax")
        choice = request(
            terminal,
            GREEN_FG("Do you want to parse the final or all atomic coordinates?") |> string,
            RadioMenu(["final", "all"]),
        )
        if choice == "final"
            ap = tryparsefinal(AtomicPositionsCard, str)
        else
            ap = tryparseall(AtomicPositionsCard, str)
        end
        if calculation == "vc-relax"
            choice = request(
                terminal,
                GREEN_FG("Do you want to parse the final or all cell parameters?") |> string,
                RadioMenu(["final", "all"]),
            )
            if choice == "final"
                cp = tryparsefinal(AtomicPositionsCard, str)
            else
                cp = tryparseall(AtomicPositionsCard, str)
            end
        end
    end
end # function output_parser

end
