module PHonon

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO.Namelists.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist
using QuantumESPRESSO.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput

using ...Namelists: namelist_helper
using ..Inputs

function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:PhInput}
    ph = namelist_helper(terminal, PhNamelist)
end
function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:Q2rInput}
    return Q2rInput(namelist_helper(terminal, Q2rNamelist))
end
function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:MatdynInput}
    ph = namelist_helper(terminal, MatdynNamelist)
end
function Inputs.input_helper(terminal::TTYTerminal, ::Type{T}) where {T<:DynmatInput}
    return DynmatInput(namelist_helper(terminal, DynmatNamelist))
end

end # module PHonon
