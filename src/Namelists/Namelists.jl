module Namelists

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using QuantumESPRESSO: to_qe
using QuantumESPRESSO.Namelists: Namelist
using Setfield: PropertyLens, set

using ..Wizard: color_string, @c_str

export namelist_helper

function namelist_helper end

# This is a helper function and should not be exported.
function setfield_helper(terminal::TTYTerminal, nml::T) where {T<:Namelist}
    while true
        print(terminal, color_string(to_qe(nml), 'b'))
        # It will continuously print until the user chooses `"no"`, i.e., he/she is satisfied.
        isdone = pairs((false, true))[request(
            terminal,
            color_string(
                "We generate an example `$(nameof(T))`. Want to change/add any field?",
                'r',
            ),
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            while true
                print(terminal, c"Type a field name: "r)
                field = strip(readline(terminal)) |> Symbol
                # Once a field successfully changes, go back to the above menu.
                # The code will asks the user whether to change another field.
                if hasfield(T, field)
                    print(terminal, c"Type its value: "r)
                    try
                        nml = set(
                            nml,
                            PropertyLens{field}(),
                            parse(fieldtype(T, field), readline(terminal)),
                        )
                    catch e
                        !isa(e, AssertionError) && rethrow(e)
                        println(terminal, c"A wrong value is given to! Try a new one!"g)
                        continue
                    end
                    break
                end
                # If the field has a wrong name, go back to `"Type a field name: "`.
                println(terminal, c"Unknown field given! Try again!"r)
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
