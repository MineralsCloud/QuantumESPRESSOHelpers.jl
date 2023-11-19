using AbInitioSoftwareBase: Namelist
using Accessors: PropertyLens, set
using REPL.TerminalMenus: request, terminal
using Term: @green, @red

export InputBuilder, build

struct InputBuilder <: Helper end

(builder::InputBuilder)(nml::Type{<:Namelist}) = builder(terminal, typeof(nml))

const build = InputBuilder()

struct FieldSetter <: Helper end

(setter::FieldSetter)(nml::Type{<:Namelist}) = setter(terminal, typeof(nml))

# This is a helper function and should not be exported.
function (setter::FieldSetter)(io::IO, nml::Namelist)
    while true
        user_response = request(io, @green("Want to change/add any field?"), YES_NO_MENU)
        println(io, @green "I'll print the result for you:")
        println(io, @green string(nml))
        if Bool(only(user_response) - 1)
            break  # Break the loop if user chooses 'no'
        end
        nml = _set(io, nml)
    end
    return nml
end

function _set(io, nml)
    fields = string.(fieldnames(typeof(nml)))
    try
        println(io, @green "Allowed fields: $fields")
        print(io, @green "Type a field name: ")
        field = Symbol(strip(readline(io)))
        if hasfield(typeof(nml), field)
            print(io, @green "Type its value: ")
            S = fieldtype(typeof(nml), field)
            nml = if S <: AbstractString
                set(nml, PropertyLens{field}(), chomp(readline(io)))
            else
                set(nml, PropertyLens{field}(), parse(S, readline(io)))
            end
        else
            println(io, @red "Unknown field given! Try again!")
        end
    catch e
        if e isa AssertionError
            println(io, @red "A wrong value is given! Try a new one!")
        elseif e isa InterruptException
        else
            rethrow(e)
        end
    end
    return nml
end
