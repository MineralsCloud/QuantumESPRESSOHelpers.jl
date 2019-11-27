module Namelists

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: BLUE_FG, GREEN_FG, RED_FG
using QuantumESPRESSOBase: to_qe
using QuantumESPRESSOBase.Namelists: Namelist
using Setfield: PropertyLens, set

export namelist_helper

function namelist_helper end

# This is a helper function and should not be exported.
function setfield_helper(terminal::TTYTerminal, nml::T) where {T<:Namelist}
    while true
        print(terminal, BLUE_FG(to_qe(nml)) |> string)
        # It will continuously print until the user chooses `"no"`, i.e., he/she is satisfied.
        isdone = pairs((false, true))[request(
            terminal,
            GREEN_FG("We generate an example `$(nameof(T))`. Want to change/add any field?") |> string,
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            while true
                print(terminal, GREEN_FG("Type a field name: ") |> string)
                field = strip(readline(terminal)) |> Symbol
                # Once a field successfully changes, go back to the above menu.
                # The code will asks the user whether to change another field.
                if hasfield(T, field)
                    print(terminal, GREEN_FG("Type its value: ") |> string)
                    try
                        nml = set(
                            nml,
                            PropertyLens{field}(),
                            parse(fieldtype(T, field), readline(terminal)),
                        )
                    catch e
                        !isa(e, AssertionError) && rethrow(e)
                        println(terminal, RED_FG("A wrong value is given to! Try a new one!") |> string)
                        continue
                    end
                    break
                end
                # If the field has a wrong name, go back to `"Type a field name: "`.
                println(terminal, RED_FG("Unknown field given! Try again!") |> string)
                continue
            end
            continue
        end
        break
    end
    return nml
end # function change_field_helper

include("PWscf.jl")
include("PHonon.jl")

end
