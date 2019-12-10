module PHonon

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSOBase.Namelists.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist
using QuantumESPRESSOBase.Inputs.PHonon: PhInput, Q2rInput, MatdynInput, DynmatInput

using ...Namelists: namelist_builder
using ..Inputs

function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:PhInput}
    ph = namelist_builder(terminal, PhNamelist)
end
function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:Q2rInput}
    return Q2rInput(namelist_builder(terminal, Q2rNamelist))
end
function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:MatdynInput}
    ph = namelist_builder(terminal, MatdynNamelist)
end
function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:DynmatInput}
    return DynmatInput(namelist_builder(terminal, DynmatNamelist))
end

end # module PHonon