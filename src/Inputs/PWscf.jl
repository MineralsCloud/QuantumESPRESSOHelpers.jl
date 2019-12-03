module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: RED_FG
using QuantumESPRESSOBase: asfieldname
using QuantumESPRESSOBase.Namelists.PWscf:
    ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist
using QuantumESPRESSOBase.Cards.PWscf: AtomicSpecies, AtomicSpeciesCard, AtomicPosition, AtomicPositionsCard, KPointsCard, CellParametersCard
using QuantumESPRESSOBase.Inputs.PWscf: PWInput

using ...Namelists: namelist_helper
using ...Cards: card_helper
using ..Inputs

function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:PWInput}
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, asfieldname(S) => namelist_helper(terminal, S))
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
                push!(fields, asfieldname(IonsNamelist) => namelist_helper(terminal, IonsNamelist))
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
                push!(fields, asfieldname(CellNamelist) => namelist_helper(terminal, CellNamelist))
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
            push!(fields, asfieldname(KPointsCard) => card_helper(terminal, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(terminal, RED_FG("Something wrong happens, try again!") |> string)
        end
    end
    push!(fields, asfieldname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(fields, asfieldname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]))
    return T(; fields...)
end # function input_builder

end # module PWscf
