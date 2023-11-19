using DataFrames: DataFrame
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

using ..QuantumESPRESSOHelpers: YES_NO_MENU, QuantumESPRESSOHelper

export OutputParser, parse_output

const FINAL_ALL_MENU = RadioMenu(Base.vect("final", "all"); charset=:ascii)

struct OutputParser <: QuantumESPRESSOHelper end

function (::OutputParser)(term::IO)
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
    calculation = CALCULATION[request(
        term,
        @green("What exact calculation is this output?"),
        RadioMenu(CALCULATION; charset=:ascii),
    )]
    if calculation == "relax" || calculation == "vc-relax"
        choice = request(
            term,
            @green("Do you want to parse the final or all atomic positions?"),
            FINAL_ALL_MENU,
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
                FINAL_ALL_MENU,
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
    choice = request(term, @green("Do you want to parse its summary?"), YES_NO_MENU)
    if choice == 1
        preamble = parse(Preamble, str)
        println(term, DataFrame(Base.vect(preamble)))
    end
    choice = request(term, @green("Do you want to parse the energies?"), YES_NO_MENU)
    if choice == 1
        df = DataFrame(eachconvergedenergy(str))
        println(term, df)
    end
    choice = request(term, @green("Do you want to parse the time used?"), YES_NO_MENU)
    if choice == 1
        df = DataFrame(eachtimeditem(str))
        println(term, df)
    end
    choice = request(
        term, @green("Do you want to parse the diagonalization info?"), YES_NO_MENU
    )
    if choice == 1
        df = DataFrame(eachdiagonalization(str))
        println(term, df)
    end
end

const parse_output = OutputParser()
