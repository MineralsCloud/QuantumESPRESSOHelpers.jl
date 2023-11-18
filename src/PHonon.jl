module PHonon

using Crayons.Box: GREEN_FG
using QuantumESPRESSOBase.PHonon:
    PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist, Q2rInput, DynmatInput
using REPL.TerminalMenus: RadioMenu, request

function build(term::IO, ::Type{PhNamelist})
    print(
        term,
        string(
            GREEN_FG(
                "Please input the atomic mass [amu] of each atomic type `amass` (separated by spaces): ",
            ),
        ),
    )
    amass = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    epsil_pool = pairs((false, true))
    epsil = epsil_pool[request(
        term, string(GREEN_FG("Please select the `epsil`: ")), RadioMenu([false, true])
    )]
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        term,
        string(GREEN_FG("Please select the `q_in_band_form`: ")),
        RadioMenu([false, true]),
    )]
    print(
        term,
        string(
            GREEN_FG(
                "Please input parameters of the Monkhorst-Pack grid `nq` 1-3 (separated by spaces): ",
            ),
        ),
    )
    nq1, nq2, nq3 = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    print(
        term,
        string(
            GREEN_FG(
                "Please input parameters of the Monkhorst-Pack grid `nk` 1-3 (separated by spaces): ",
            ),
        ),
    )
    nk1, nk2, nk3 = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    print(term, string(GREEN_FG("Please input offset `k` 1-3 (separated by spaces): ")))
    k1, k2, k3 = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
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
    return help_set(term, ph)
end
function build(term::IO, ::Type{T}) where {T<:Q2rNamelist}
    print(term, string(GREEN_FG("name of input dynamical matrices `fildyn`: ")))
    fildyn = strip(readline(term))
    print(term, string(GREEN_FG("name of output force constants `flfrc`: ")))
    flfrc = strip(readline(term))
    zasr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    zasr = zasr_pool[request(
        term,
        string(
            GREEN_FG(
                "Please input the type of acoustic sum rules used for the Born effective charges `zasr`: ",
            ),
        ),
        RadioMenu(collect(values(zasr_pool))),
    )]
    q2r = T(; fildyn=fildyn, flfrc=flfrc, zasr=zasr)
    return help_set(term, q2r)
end
function build(term::IO, ::Type{T}) where {T<:MatdynNamelist}
    dos_pool = pairs((false, true))
    dos = dos_pool[request(
        term,
        string(GREEN_FG("Please select if calculate phonon density of states `dos`: ")),
        RadioMenu([false, true]),
    )]
    print(term, string(GREEN_FG("Please input the energy step, in cm^(-1) `deltaE`: ")))
    deltaE = parse(Float64, readline(term))
    print(
        term,
        string(
            GREEN_FG(
                "Please input uniform q-point grid for DOS calculation `nk` 1-3 (separated by spaces): ",
            ),
        ),
    )
    nk1, nk2, nk3 = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        term,
        string(GREEN_FG("Please input the type of acoustic sum rule `asr`: ")),
        RadioMenu(collect(values(asr_pool))),
    )]
    print(term, string(GREEN_FG("name of output force constants `flfrc`: ")))
    flfrc = strip(readline(term))
    print(term, string(GREEN_FG("name of input dynamical matrices `fildyn`: ")))
    fildyn = strip(readline(term))
    print(
        term,
        string(
            GREEN_FG(
                "Please input the masses of atoms in the supercell (a.m.u.) `amass` (separated by spaces): ",
            ),
        ),
    )
    amass = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    print(
        term,
        string(GREEN_FG("Please input the number of atom types in the supercell `ntyp`: ")),
    )
    ntyp = parse(Int, readline(term))
    q_in_band_form_pool = pairs((false, true))
    q_in_band_form = q_in_band_form_pool[request(
        term,
        string(GREEN_FG("Please select the `q_in_band_form`: ")),
        RadioMenu([false, true]),
    )]
    q_in_cryst_coord_pool = pairs((false, true))
    q_in_cryst_coord = q_in_cryst_coord_pool[request(
        term,
        string(GREEN_FG("Please select the `q_in_cryst_coord`: ")),
        RadioMenu([false, true]),
    )]
    nosym_pool = pairs((false, true))
    nosym = nosym_pool[request(
        term,
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
    return help_set(term, matdyn)
end
function build(term::IO, ::Type{T}) where {T<:DynmatNamelist}
    asr_pool = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
    asr = asr_pool[request(
        term,
        string(GREEN_FG("Please select the type of acoustic sum rule `asr`: ")),
        RadioMenu(collect(values(asr_pool))),
    )]
    print(
        term,
        string(
            GREEN_FG("Please input mass for each atom type `amass` (separated by spaces): ")
        ),
    )
    amass = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    dynmat = T(; asr=asr, amass=amass)
    return help_set(term, dynmat)
end

# function build(term::IO, ::Type{PhInput})
#     return build(term, PhNamelist)
# end
function build(term::IO, ::Type{Q2rInput})
    return Q2rInput(build(term, Q2rNamelist))
end
# function build(term::IO, ::Type{MatdynInput})
#     return build(term, MatdynNamelist)
# end
function build(term::IO, ::Type{DynmatInput})
    return DynmatInput(build(term, DynmatNamelist))
end

end
