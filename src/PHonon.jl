module PHonon

using QuantumESPRESSOBase.PHonon:
    PhNamelist,
    Q2rNamelist,
    MatdynNamelist,
    DynmatNamelist,
    PhInput,
    Q2rInput,
    MatdynInput,
    DynmatInput
using REPL.TerminalMenus: RadioMenu, request
using Term: @green

using ..QuantumESPRESSOHelpers: YES_NO_MENU, help_set

const EPSIL = pairs((false, true))
const Q_IN_BAND_FORM = pairs((false, true))
const ZASR = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
const DOS = pairs((false, true))
const ASR = pairs(("no", "simple", "crystal", "one-dim", "zero-dim"))
const Q_IN_CRYST_COORD = pairs((false, true))
const NOSYM = pairs((false, true))

function build(term::IO, ::Type{PhNamelist})
    print(
        term,
        @green "Please input the atomic mass [amu] of each atomic type `amass` (separated by spaces): "
    )
    amass = map(Base.Fix1(parse, Float64), split(readline(term), " "; keepempty=false))
    epsil = EPSIL[request(term, @green("Please select the `epsil`: "), YES_NO_MENU)]
    q_in_band_form = Q_IN_BAND_FORM[request(
        term, @green("Please select the `q_in_band_form`: "), YES_NO_MENU
    )]
    print(
        term,
        @green "Please input parameters of the Monkhorst-Pack grid `nq` 1-3 (separated by spaces): "
    )
    nq1, nq2, nq3 = map(
        Base.Fix1(parse, Float64), split(readline(term), " "; keepempty=false)
    )
    print(
        term,
        @green "Please input parameters of the Monkhorst-Pack grid `nk` 1-3 (separated by spaces): "
    )
    nk1, nk2, nk3 = map(Base.Fix1(parse, Int), split(readline(term), " "; keepempty=false))
    print(term, @green "Please input offset `k` 1-3 (separated by spaces): ")
    k1, k2, k3 = map(Base.Fix1(parse, Int), split(readline(term), " "; keepempty=false))
    ph = PhNamelist(;
        amass, epsil, q_in_band_form, nq1, nq2, nq3, nk1, nk2, nk3, k1, k2, k3
    )
    return help_set(term, ph)
end
function build(term::IO, ::Type{Q2rNamelist})
    print(term, @green "name of input dynamical matrices `fildyn`: ")
    fildyn = strip(readline(term))
    print(term, @green "name of output force constants `flfrc`: ")
    flfrc = strip(readline(term))
    zasr = ZASR[request(
        term,
        @green(
            "Please input the type of acoustic sum rules used for the Born effective charges `zasr`: "
        ),
        RadioMenu(ZASR; charset=:ascii),
    )]
    q2r = Q2rNamelist(; fildyn, flfrc, zasr)
    return help_set(term, q2r)
end
function build(term::IO, ::Type{MatdynNamelist})
    dos = DOS[request(
        term,
        @green("Please select if calculate phonon density of states `dos`: "),
        YES_NO_MENU,
    )]
    print(term, @green "Please input the energy step, in cm^(-1) `deltaE`: ")
    deltaE = parse(Float64, readline(term))
    print(
        term,
        @green "Please input uniform q-point grid for DOS calculation `nk` 1-3 (separated by spaces): "
    )
    nk1, nk2, nk3 = map(Base.Fix1(parse, Int), split(readline(term), " "; keepempty=false))
    asr = ASR[request(
        term,
        @green("Please input the type of acoustic sum rule `asr`: "),
        RadioMenu(ASR; charset=:ascii),
    )]
    print(term, @green "name of output force constants `flfrc`: ")
    flfrc = strip(readline(term))
    print(term, @green "name of input dynamical matrices `fildyn`: ")
    fildyn = strip(readline(term))
    print(
        term,
        @green "Please input the masses of atoms in the supercell (a.m.u.) `amass` (separated by spaces): "
    )
    amass = map(Base.Fix1(parse, Float64), split(readline(term), " "; keepempty=false))
    print(term, @green "Please input the number of atom types in the supercell `ntyp`: ")
    ntyp = parse(Int, readline(term))
    q_in_band_form = Q_IN_BAND_FORM[request(
        term, @green("Please select the `q_in_band_form`: "), YES_NO_MENU
    )]
    q_in_cryst_coord = Q_IN_CRYST_COORD[request(
        term, @green("Please select the `q_in_cryst_coord`: "), YES_NO_MENU
    )]
    nosym = NOSYM[request(
        term,
        @green("Please select if impose symmetry and time reversal `nosym`: "),
        YES_NO_MENU,
    )]
    matdyn = MatdynNamelist(;
        dos,
        deltaE,
        nk1,
        nk2,
        nk3,
        asr,
        flfrc,
        fildyn,
        amass,
        ntyp,
        q_in_band_form,
        q_in_cryst_coord,
        nosym,
    )
    return help_set(term, matdyn)
end
function build(term::IO, ::Type{DynmatNamelist})
    asr = ASR[request(
        term,
        @green("Please select the type of acoustic sum rule `asr`: "),
        RadioMenu(ASR; charset=:ascii),
    )]
    print(
        term, @green "Please input mass for each atom type `amass` (separated by spaces): "
    )
    amass = map(Base.Fix1(parse, Float64), split(readline(term), " "; keepempty=false))
    dynmat = DynmatNamelist(; asr, amass)
    return help_set(term, dynmat)
end
build(term::IO, ::Type{PhInput}) = PhInput(build(term, PhNamelist))
build(term::IO, ::Type{Q2rInput}) = Q2rInput(build(term, Q2rNamelist))
build(term::IO, ::Type{MatdynInput}) = MatdynInput(build(term, MatdynNamelist))
build(term::IO, ::Type{DynmatInput}) = DynmatInput(build(term, DynmatNamelist))

end
