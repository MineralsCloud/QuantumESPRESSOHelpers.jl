using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: BLUE_FG, GREEN_FG, RED_FG
using QuantumESPRESSOBase.Inputs: Namelist, qestring
using Setfield: PropertyLens, set

export build

function build end

# This is a helper function and should not be exported.
function setfield_helper(terminal::TTYTerminal, nml::T) where {T<:Namelist}
    while true
        print(terminal, string(BLUE_FG(qestring(nml))))
        # It will continuously print until the user chooses `"no"`, i.e., he/she is satisfied.
        isdone = pairs((false, true))[request(
            terminal,
            string(
                GREEN_FG(
                    "We generate an example `$(nameof(T))`. Want to change/add any field?"
                ),
            ),
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            while true
                print(terminal, string(GREEN_FG("Type a field name: ")))
                field = Symbol(strip(readline(terminal)))
                # Once a field successfully changes, go back to the above menu.
                # The code will asks the user whether to change another field.
                if hasfield(T, field)
                    print(terminal, string(GREEN_FG("Type its value: ")))
                    try
                        S = fieldtype(T, field)
                        nml = if S <: AbstractString
                            set(nml, PropertyLens{field}(), chomp(readline(terminal)))
                        else
                            set(
                                nml,
                                PropertyLens{field}(),
                                parse(fieldtype(T, field), readline(terminal)),
                            )
                        end
                    catch e
                        !isa(e, AssertionError) && rethrow(e)
                        println(
                            terminal,
                            string(RED_FG("A wrong value is given to! Try a new one!")),
                        )
                        continue
                    end
                    break
                end
                # If the field has a wrong name, go back to `"Type a field name: "`.
                println(terminal, string(RED_FG("Unknown field given! Try again!")))
                continue
            end
            continue
        end
        break
    end
    return nml
end

include("PWscf.jl")
include("PHonon.jl")
