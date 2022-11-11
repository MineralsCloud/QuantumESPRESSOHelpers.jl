module PHonon

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: GREEN_FG
using QuantumESPRESSOBase.Inputs.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist, Q2rInput, DynmatInput

using ..Inputs: setfield_helper

import ..Inputs: namelist_builder, input_builder

function namelist_builder(terminal::TTYTerminal, ::Type{T}) where {T<:PhNamelist}
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the atomic mass [amu] of each atomic type `amass` (separated by spaces): ",
            ),
        ),
    )
    amass = map(x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false))
    epsil_pool = pairs((false, true))
    epsil = epsil_pool[request(
        terminal, string(GREEN_FG("Please select the `epsil`: ")), RadioMenu([false, true])
    )]
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        terminal,
        string(GREEN_FG("Please select the `q_in_band_form`: ")),
        RadioMenu([false, true]),
    )]
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input parameters of the Monkhorst-Pack grid `nq` 1-3 (separated by spaces): ",
            ),
        ),
    )
    nq1, nq2, nq3 = map(
        x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false)
    )
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input parameters of the Monkhorst-Pack grid `nk` 1-3 (separated by spaces): ",
            ),
        ),
    )
    nk1, nk2, nk3 = map(
        x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false)
    )
    print(terminal, string(GREEN_FG("Please input offset `k` 1-3 (separated by spaces): ")))
    k1, k2, k3 = map(
        x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false)
    )
    ph = T(;
        amass=amass,
        epsil=epsil,
        q_in_band_form=q_in_band_form,
        nq1=nq1,
        nq2=nq2,
        nq3=nq3,
        nk1=nk1,
        nk2=nk2,
        nk3=nk3,
        k1=k1,
        k2=k2,
        k3=k3,
    )
    return setfield_helper(terminal, ph)
end
function namelist_builder(terminal::TTYTerminal, ::Type{T}) where {T<:Q2rNamelist}
    print(terminal, string(GREEN_FG("name of input dynamical matrices `fildyn`: ")))
    fildyn = strip(readline(terminal))
    print(terminal, string(GREEN_FG("name of output force constants `flfrc`: ")))
    flfrc = strip(readline(terminal))
    zasr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    zasr = zasr_pool[request(
        terminal,
        string(
            GREEN_FG(
                "Please input the type of acoustic sum rules used for the Born effective charges `zasr`: ",
            ),
        ),
        RadioMenu(collect(values(zasr_pool))),
    )]
    q2r = T(; fildyn=fildyn, flfrc=flfrc, zasr=zasr)
    return setfield_helper(terminal, q2r)
end
function namelist_builder(terminal::TTYTerminal, ::Type{T}) where {T<:MatdynNamelist}
    dos_pool = pairs((false, true))
    dos = dos_pool[request(
        terminal,
        string(GREEN_FG("Please select if calculate phonon density of states `dos`: ")),
        RadioMenu([false, true]),
    )]
    print(terminal, string(GREEN_FG("Please input the energy step, in cm^(-1) `deltaE`: ")))
    deltaE = parse(Float64, readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input uniform q-point grid for DOS calculation `nk` 1-3 (separated by spaces): ",
            ),
        ),
    )
    nk1, nk2, nk3 = map(
        x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false)
    )
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        terminal,
        string(GREEN_FG("Please input the type of acoustic sum rule `asr`: ")),
        RadioMenu(collect(values(asr_pool))),
    )]
    print(terminal, string(GREEN_FG("name of output force constants `flfrc`: ")))
    flfrc = strip(readline(terminal))
    print(terminal, string(GREEN_FG("name of input dynamical matrices `fildyn`: ")))
    fildyn = strip(readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the masses of atoms in the supercell (a.m.u.) `amass` (separated by spaces): ",
            ),
        ),
    )
    amass = map(x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false))
    print(
        terminal,
        string(GREEN_FG("Please input the number of atom types in the supercell `ntyp`: ")),
    )
    ntyp = parse(Int, readline(terminal))
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        terminal,
        string(GREEN_FG("Please select the `q_in_band_form`: ")),
        RadioMenu([false, true]),
    )]
    q_in_cryst_coord_pool = pairs((false, true))
    q_in_cryst_coord = q_in_cryst_coord_pool[request(
        terminal,
        string(GREEN_FG("Please select the `q_in_cryst_coord`: ")),
        RadioMenu([false, true]),
    )]
    nosym_pool = pairs((false, true))
    nosym = nosym_pool[request(
        terminal,
        string(GREEN_FG("Please select if impose symmetry and time reversal `nosym`: ")),
        RadioMenu([false, true]),
    )]
    matdyn = T(;
        dos=dos,
        deltaE=deltaE,
        nk1=nk1,
        nk2=nk2,
        nk3=nk3,
        asr=asr,
        flfrc=flfrc,
        fildyn=fildyn,
        amass=amass,
        ntyp=ntyp,
        q_in_band_form=q_in_band_form,
        q_in_cryst_coord=q_in_cryst_coord,
        nosym=nosym,
    )
    return setfield_helper(terminal, matdyn)
end
function namelist_builder(terminal::TTYTerminal, ::Type{T}) where {T<:DynmatNamelist}
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        terminal,
        string(GREEN_FG("Please select the type of acoustic sum rule `asr`: ")),
        RadioMenu(collect(values(asr_pool))),
    )]
    print(
        terminal,
        string(
            GREEN_FG("Please input mass for each atom type `amass` (separated by spaces): ")
        ),
    )
    amass = map(x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false))
    dynmat = T(; asr=asr, amass=amass)
    return setfield_helper(terminal, dynmat)
end

# function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:PhInput}
#     ph = namelist_builder(terminal, PhNamelist)
# end
function input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:Q2rInput}
    return Q2rInput(namelist_builder(terminal, Q2rNamelist))
end
# function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:MatdynInput}
#     ph = namelist_builder(terminal, MatdynNamelist)
# end
function input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:DynmatInput}
    return DynmatInput(namelist_builder(terminal, DynmatNamelist))
end

end
