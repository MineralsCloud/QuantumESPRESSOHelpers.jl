using Parameters: type2dict
using PrettyTables: Highlighter, pretty_table, ft_printf
using QuantumESPRESSOBase.PWscf: AtomicPositionsCard, CellParametersCard
using QuantumESPRESSOParser.PWscf
using Term: @blue

export PWOutput, output_parser

struct PWOutput end

function output_parser(term::IO, ::Type{T}) where {T<:PWOutput}
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    print(term, @green "Please give the absolute path to your output file: ")
    path = abspath(strip(readline(term)))
    str = try
        open(path, "r") do io
            read(io, String)
        end
    catch
        @warn("File '$path' not found!")
    end
    if isjobdone(str)
        println(term, @green "The job is done! Ready to parse!")
    else
        println(term, @red "The job is not finished, be careful!")
    end
    calculation = calculations[request(
        term,
        @green("What exact calculation is this output?"),
        RadioMenu(collect(values(calculations))),
    )]
    if calculation ∈ ("relax", "vc-relax")
        choice = request(
            term,
            @green("Do you want to parse the final or all atomic positions?"),
            RadioMenu(["final", "all"]),
        )
        if choice == 1
            ap = tryparsefinal(AtomicPositionsCard, str)
            println(term, @blue "Print the final atomic positions:")
            display(ap.data)
        else
            ap = tryparseall(AtomicPositionsCard, str)
            println(term, @blue "Print all atomic positions:")
            foreach(x -> display(x.data), ap)
        end
        if calculation == "vc-relax"
            choice = request(
                term,
                @green("Do you want to parse the final or all cell parameters?"),
                RadioMenu(["final", "all"]),
            )
            if choice == 1
                cp = tryparsefinal(CellParametersCard, str)
                println(term, @blue "Print the final cell parameters:")
                println(term, cp.data)
            else
                cp = tryparseall(CellParametersCard, str)
                println(term, @blue "Print all cell parameters:")
                foreach(x -> display(x.data), cp)
            end
            if isrelaxed(str)
                println(term, @green "The structure is relaxed!")
            else
                println(term, @red "The structure is not well-relaxed!")
            end
        end
    end
    hl_odd = Highlighter(; f=(data, i, j) -> (i % 2) == 0)
    choice = request(
        term, @green("Do you want to parse its summary?"), RadioMenu(["yes", "no"])
    )
    if choice == 1
        preamble = parse(Preamble, str)
        pretty_table(type2dict(preamble); highlighters=hl_odd)
    end
    choice = request(
        term, @green("Do you want to parse the energies?"), RadioMenu(["yes", "no"])
    )
    if choice == 1
        df = parse_electrons_energies(str, :combined)
        pretty_table(df; highlighters=hl_odd, formatter=ft_printf("%10.5f"))
    end
    choice = request(
        term, @green("Do you want to parse the time used?"), RadioMenu(["yes", "no"])
    )
    if choice == 1
        df = parse_clock(str)
        pretty_table(df; highlighters=hl_odd, formatter=ft_printf("%10.5f"))
    end
    choice = request(
        term,
        @green("Do you want to parse the diagonalization info?"),
        RadioMenu(["yes", "no"]),
    )
    if choice == 1
        df = Outputs.PWscf.parse_diagonalization(str)
        pretty_table(df; highlighters=hl_odd, formatter=ft_printf("%10.5f"))
    end
end
