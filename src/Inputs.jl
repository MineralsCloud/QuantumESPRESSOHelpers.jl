using AbInitioSoftwareBase: Namelist
using Accessors: PropertyLens, set
using Crayons.Box: BLUE_FG
using REPL.TerminalMenus: RadioMenu, request, terminal
using Term: @green, @red

export build

function build end

# This is a helper function and should not be exported.
function help_set(term, nml::Namelist)
    while true
        print(term, @green string(nml))
        user_response = request(
            term,
            @green(
                "the current namelist is `$(nameof(typeof(nml)))`. Want to change/add any field?"
            ),
            RadioMenu(Base.vect("yes", "no"); charset=:ascii),
        )
        if only(user_response) - 1
            break  # Break the loop if user chooses 'no'
        end
        nml = _help_set_iter(term, nml)
    end
    return nml
end
help_set(nml) = help_set(terminal, nml)

function _help_set_iter(term, nml::Namelist)
    try
        while true
            print(term, @green "Type a field name: ")
            field = Symbol(strip(readline(term)))
            if hasfield(typeof(nml), field)
                print(term, @green "Type its value: ")
                try
                    S = fieldtype(typeof(nml), field)
                    nml = if S <: AbstractString
                        set(nml, PropertyLens{field}(), chomp(readline(term)))
                    else
                        set(nml, PropertyLens{field}(), parse(S, readline(term)))
                    end
                catch e
                    if !(e isa AssertionError)
                        rethrow(e)
                    end
                    println(term, @red "A wrong value is given! Try a new one!")
                end
            else
                println(term, @red "Unknown field given! Try again!")
            end
        end
    catch e
        if !(e isa InterruptException)
            rethrow(e)
        end
    end
    return nml
end

include("PWscf.jl")
include("PHonon.jl")
