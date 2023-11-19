module QuantumESPRESSOHelpers

using REPL.TerminalMenus: RadioMenu

const YES_NO_MENU = RadioMenu(Base.vect("yes", "no"); charset=:ascii)

abstract type Helper end

include("input.jl")
include("PWscf/PWscf.jl")
include("PHonon.jl")

end
