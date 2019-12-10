module Cards

export card_builder

function card_builder end

module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: GREEN_FG
using QuantumESPRESSOBase.Cards.PWscf: GammaPoint, MonkhorstPackGrid, KPointsCard

using ..Cards

function Cards.card_builder(terminal::TTYTerminal, ::Type{T}) where {T<:KPointsCard}
    kpt_style = request(
        terminal,
        GREEN_FG("What k-point style do you want?") |> string,
        RadioMenu(["gamma", "automatic"]),
    )
    return if kpt_style == 1
        KPointsCard("gamma", GammaPoint())
    else  # "automatic"
        print(terminal, GREEN_FG("What 3-element k-point grid do you want (separated by spaces): ") |> string)
        grid = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        print(terminal, GREEN_FG("What 3-element k-point offsets do you want (separated by spaces): ") |> string)
        offsets = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        return KPointsCard("automatic", MonkhorstPackGrid(grid, offsets))
    end
end # function card_builder

end # module PWscf

end
