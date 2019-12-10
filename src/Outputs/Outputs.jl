module Outputs

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons: Crayon
using Crayons.Box: GREEN_FG, BLUE_FG
using Parameters: type2dict
using PrettyTables: Highlighter, pretty_table, ft_printf
using QuantumESPRESSOBase.Cards.PWscf: AtomicPositionsCard, CellParametersCard
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
            GREEN_FG("Do you want to parse the final or all atomic positions?") |> string,
            RadioMenu(["final", "all"]),
        )
        if choice == 1
            ap = tryparsefinal(AtomicPositionsCard, str)
            println(terminal, BLUE_FG("Print the final atomic positions:"))
            display(ap.data)
        else
            ap = tryparseall(AtomicPositionsCard, str)
            println(terminal, BLUE_FG("Print all atomic positions:"))
            foreach(x -> display(x.data), ap)
        end
        if calculation == "vc-relax"
            choice = request(
                terminal,
                GREEN_FG("Do you want to parse the final or all cell parameters?") |> string,
                RadioMenu(["final", "all"]),
            )
            if choice == 1
                cp = tryparsefinal(CellParametersCard, str)
                println(terminal, BLUE_FG("Print the final cell parameters:"))
                println(terminal, cp.data)
            else
                cp = tryparseall(CellParametersCard, str)
                println(terminal, BLUE_FG("Print all cell parameters:"))
                foreach(x -> display(x.data), cp)
            end
        end
    end
    hl_odd = Highlighter(
        f = (data, i, j) -> (i % 2) == 0,
        crayon = Crayon(background = :blue)
    )
    choice = request(
        terminal,
        GREEN_FG("Do you want to parse its summary?") |> string,
        RadioMenu(["yes", "no"]),
    )
    if choice == 1
        preamble = parse(Preamble, str)
        pretty_table(preamble |> type2dict; highlighters = hl_odd)
    end
    choice = request(
        terminal,
        GREEN_FG("Do you want to parse the energies?") |> string,
        RadioMenu(["yes", "no"]),
    )
    if choice == 1
        df = parse_electrons_energies(str, :combined)
        pretty_table(df; highlighters = hl_odd, formatter = ft_printf("%10.5f"))
    end
    choice = request(
        terminal,
        GREEN_FG("Do you want to parse the time used?") |> string,
        RadioMenu(["yes", "no"]),
    )
    if choice == 1
        df = parse_clock(str)
        pretty_table(df; highlighters = hl_odd, formatter = ft_printf("%10.5f"))
    end
    choice = request(
        terminal,
        GREEN_FG("Do you want to parse the diagonalization info?") |> string,
        RadioMenu(["yes", "no"]),
    )
    if choice == 1
        df = Outputs.PWscf.parse_diagonalization(str)
        pretty_table(df; highlighters = hl_odd, formatter = ft_printf("%10.5f"))
    end
end # function output_parser

end
