using AbInitioSoftwareBase: Namelist
using Crayons.Box: BLUE_FG
using REPL.TerminalMenus: RadioMenu, request, terminal
using Setfield: PropertyLens, set
using Term: @green, @red

export build

function build end

# This is a helper function and should not be exported.
function setfield_helper(term::IO, nml::T) where {T<:Namelist}
    while true
        print(term, @green string(nml))
        # It will continuously print until the user chooses `"no"`, i.e., he/she is satisfied.
        isdone = Base.vect(false, true)[request(
            term,
            @green("We generate an example `$(nameof(T))`. Want to change/add any field?"),
            RadioMenu(["yes", "no"]),
        )]
        if !isdone
            while true
                print(term, @green "Type a field name: ")
                field = Symbol(strip(readline(term)))
                # Once a field successfully changes, go back to the above menu.
                # The code will asks the user whether to change another field.
                if hasfield(T, field)
                    print(term, @green "Type its value: ")
                    try
                        S = fieldtype(T, field)
                        nml = if S <: AbstractString
                            set(nml, PropertyLens{field}(), chomp(readline(term)))
                        else
                            set(
                                nml,
                                PropertyLens{field}(),
                                parse(fieldtype(T, field), readline(term)),
                            )
                        end
                    catch e
                        !isa(e, AssertionError) && rethrow(e)
                        println(term, @red "A wrong value is given to! Try a new one!")
                        continue
                    end
                    break
                end
                # If the field has a wrong name, go back to `"Type a field name: "`.
                println(term, @red "Unknown field given! Try again!")
                continue
            end
            continue
        end
        break
    end
    return nml
end
setfield_helper(nml) = setfield_helper(terminal, nml)

include("PWscf.jl")
include("PHonon.jl")
