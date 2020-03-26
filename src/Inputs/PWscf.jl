module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: RED_FG, GREEN_FG
using QuantumESPRESSOBase.Inputs: entryname, qestring
using QuantumESPRESSOBase.Inputs.PWscf:
    ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    KPointsCard,
    CellParametersCard,
    PWInput

using ...Namelists: namelist_builder
using ...Cards: card_builder
using ..Inputs

function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:PWInput}
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, entryname(S) => namelist_builder(terminal, S))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, RED_FG("Something wrong happens, try again!") |> string)
            end
        end
    end
    if fields[:control].calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(
                    fields,
                    entryname(IonsNamelist) => namelist_builder(terminal, IonsNamelist),
                )
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, RED_FG("Something wrong happens, try again!") |> string)
            end
        end
    else
        IonsNamelist()
    end
    if fields[:control].calculation ∈ ("vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(
                    fields,
                    entryname(CellNamelist) => namelist_builder(terminal, CellNamelist),
                )
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, RED_FG("Something wrong happens, try again!") |> string)
            end
        end
    else
        CellNamelist()
    end
    haserror = true
    while haserror
        try
            push!(fields, entryname(KPointsCard) => card_builder(terminal, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(terminal, RED_FG("Something wrong happens, try again!") |> string)
        end
    end
    push!(fields, entryname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(
        fields,
        entryname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]),
    )
    result = T(; fields...)
    saveresult = pairs((true, false))[request(
        terminal,
        GREEN_FG("Do you want to save the generated input to file?") |> string,
        RadioMenu(["yes", "no"]),
    )]
    if saveresult
        print(terminal, GREEN_FG("Input file name: ") |> string)
        write(readline(terminal) |> chomp, qestring(result))
    end
    return result
end # function input_builder

end # module PWscf
