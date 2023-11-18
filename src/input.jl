using AbInitioSoftwareBase: Namelist
using Accessors: PropertyLens, set
using REPL.TerminalMenus: RadioMenu, request, terminal
using Term: @green, @red

export build

function build end

# This is a helper function and should not be exported.
function help_set(term, nml::Namelist)
    while true
        user_response = request(
            term,
            @green("Want to change/add any field?"),
            RadioMenu(Base.vect("yes", "no"); charset=:ascii),
        )
        println(term, @green "I'll print the result for you:")
        println(term, @green string(nml))
        if Bool(only(user_response) - 1)
            break  # Break the loop if user chooses 'no'
        end
        nml = _help_set_iter(term, nml)
    end
    return nml
end
help_set(nml) = help_set(terminal, nml)

function _help_set_iter(term, nml::Namelist)
    fields = string.(fieldnames(typeof(nml)))
    try
        println(term, @green "Allowed fields: $fields")
        print(term, @green "Type a field name: ")
        field = Symbol(strip(readline(term)))
        if hasfield(typeof(nml), field)
            print(term, @green "Type its value: ")
            S = fieldtype(typeof(nml), field)
            nml = if S <: AbstractString
                set(nml, PropertyLens{field}(), chomp(readline(term)))
            else
                set(nml, PropertyLens{field}(), parse(S, readline(term)))
            end
        else
            println(term, @red "Unknown field given! Try again!")
        end
    catch e
        if e isa AssertionError
            println(term, @red "A wrong value is given! Try a new one!")
        elseif e isa InterruptException
        else
            rethrow(e)
        end
    end
    return nml
end