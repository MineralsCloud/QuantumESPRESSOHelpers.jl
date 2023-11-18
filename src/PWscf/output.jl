using DataFrames: DataFrame
using QuantumESPRESSOBase.PWscf: AtomicPositionsCard, CellParametersCard
using QuantumESPRESSOParser.PWscf:
    Preamble,
    isjobdone,
    isoptimized,
    eachatomicpositionscard,
    eachcellparameterscard,
    eachconvergedenergy,
    eachtimeditem,
    eachdiagonalization
using Term: @blue

export PWOutput, parse_output

struct PWOutput end

function parse_output(term::IO)
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
    calculation = CALCULATIONS[request(
        term,
        @green("What exact calculation is this output?"),
        RadioMenu(collect(values(CALCULATIONS))),
    )]
    if calculation == "relax" || calculation == "vc-relax"
        choice = request(
            term,
            @green("Do you want to parse the final or all atomic positions?"),
            RadioMenu(["final", "all"]),
        )
        cards = collect(eachatomicpositionscard(str))
        if choice == 1
            println(term, @blue "Print the final atomic positions:")
            println(term, last(cards).data)
        else
            println(term, @blue "Print all atomic positions:")
            for card in cards
                println(term, card.data)
            end
        end
        if calculation == "vc-relax"
            choice = request(
                term,
                @green("Do you want to parse the final or all cell parameters?"),
                RadioMenu(["final", "all"]),
            )
            cards = collect(eachcellparameterscard(str))
            if choice == 1
                println(term, @blue "Print the final cell parameters:")
                println(term, last(cards).data)
            else
                println(term, @blue "Print all cell parameters:")
                for card in cards
                    println(term, card.data)
                end
            end
            if isoptimized(str)
                println(term, @green "The structure is relaxed!")
            else
                println(term, @red "The structure is not well-relaxed!")
            end
        end
    end
    choice = request(
        term, @green("Do you want to parse its summary?"), RadioMenu(["yes", "no"])
    )
    if choice == 1
        preamble = parse(Preamble, str)
        println(term, DataFrame(Base.vect(preamble)))
    end
    choice = request(
        term, @green("Do you want to parse the energies?"), RadioMenu(["yes", "no"])
    )
    if choice == 1
        df = DataFrame(eachconvergedenergy(str))
        println(term, df)
    end
    choice = request(
        term, @green("Do you want to parse the time used?"), RadioMenu(["yes", "no"])
    )
    if choice == 1
        df = DataFrame(eachtimeditem(str))
        println(term, df)
    end
    choice = request(
        term,
        @green("Do you want to parse the diagonalization info?"),
        RadioMenu(["yes", "no"]),
    )
    if choice == 1
        df = DataFrame(eachdiagonalization(str))
        println(term, df)
    end
end
parse_output() = parse_output(terminal)
