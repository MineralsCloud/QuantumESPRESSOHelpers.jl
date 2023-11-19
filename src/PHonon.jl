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

using ..QuantumESPRESSOHelpers: YES_NO_MENU, InputBuilder, FieldSetter

const EPSIL = Base.vect(false, true)
const Q_IN_BAND_FORM = Base.vect(false, true)
const ZASR = Base.vect("no", "simple", "crystal", "one-dim", "zero-dim")
const DOS = Base.vect(false, true)
const ASR = Base.vect("no", "simple", "crystal", "one-dim", "zero-dim")
const Q_IN_CRYST_COORD = Base.vect(false, true)
const NOSYM = Base.vect(false, true)

function (::InputBuilder)(io::IO, ::Type{PhNamelist})
    print(
        io,
        @green "Please input the atomic mass [amu] of each atomic type `amass` (separated by spaces): "
    )
    amass = map(Base.Fix1(parse, Float64), split(readline(io), " "; keepempty=false))
    epsil = EPSIL[request(io, @green("Please select the `epsil`: "), YES_NO_MENU)]
    q_in_band_form = Q_IN_BAND_FORM[request(
        io, @green("Please select the `q_in_band_form`: "), YES_NO_MENU
    )]
    print(
        io,
        @green "Please input parameters of the Monkhorst-Pack grid `nq` 1-3 (separated by spaces): "
    )
    nq1, nq2, nq3 = map(
        Base.Fix1(parse, Float64), split(readline(io), " "; keepempty=false)
    )
    print(
        io,
        @green "Please input parameters of the Monkhorst-Pack grid `nk` 1-3 (separated by spaces): "
    )
    nk1, nk2, nk3 = map(Base.Fix1(parse, Int), split(readline(io), " "; keepempty=false))
    print(io, @green "Please input offset `k` 1-3 (separated by spaces): ")
    k1, k2, k3 = map(Base.Fix1(parse, Int), split(readline(io), " "; keepempty=false))
    ph = PhNamelist(;
        amass, epsil, q_in_band_form, nq1, nq2, nq3, nk1, nk2, nk3, k1, k2, k3
    )
    return FieldSetter()(io, ph)
end
function (::InputBuilder)(io::IO, ::Type{Q2rNamelist})
    print(io, @green "name of input dynamical matrices `fildyn`: ")
    fildyn = strip(readline(io))
    print(io, @green "name of output force constants `flfrc`: ")
    flfrc = strip(readline(io))
    zasr = ZASR[request(
        io,
        @green(
            "Please input the type of acoustic sum rules used for the Born effective charges `zasr`: "
        ),
        RadioMenu(ZASR; charset=:ascii),
    )]
    q2r = Q2rNamelist(; fildyn, flfrc, zasr)
    return FieldSetter()(io, q2r)
end
function (::InputBuilder)(io::IO, ::Type{MatdynNamelist})
    dos = DOS[request(
        io,
        @green("Please select if calculate phonon density of states `dos`: "),
        YES_NO_MENU,
    )]
    print(io, @green "Please input the energy step, in cm^(-1) `deltaE`: ")
    deltaE = parse(Float64, readline(io))
    print(
        io,
        @green "Please input uniform q-point grid for DOS calculation `nk` 1-3 (separated by spaces): "
    )
    nk1, nk2, nk3 = map(Base.Fix1(parse, Int), split(readline(io), " "; keepempty=false))
    asr = ASR[request(
        io,
        @green("Please input the type of acoustic sum rule `asr`: "),
        RadioMenu(ASR; charset=:ascii),
    )]
    print(io, @green "name of output force constants `flfrc`: ")
    flfrc = strip(readline(io))
    print(io, @green "name of input dynamical matrices `fildyn`: ")
    fildyn = strip(readline(io))
    print(
        io,
        @green "Please input the masses of atoms in the supercell (a.m.u.) `amass` (separated by spaces): "
    )
    amass = map(Base.Fix1(parse, Float64), split(readline(io), " "; keepempty=false))
    print(io, @green "Please input the number of atom types in the supercell `ntyp`: ")
    ntyp = parse(Int, readline(io))
    q_in_band_form = Q_IN_BAND_FORM[request(
        io, @green("Please select the `q_in_band_form`: "), YES_NO_MENU
    )]
    q_in_cryst_coord = Q_IN_CRYST_COORD[request(
        io, @green("Please select the `q_in_cryst_coord`: "), YES_NO_MENU
    )]
    nosym = NOSYM[request(
        io,
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
    return FieldSetter()(io, matdyn)
end
function (::InputBuilder)(io::IO, ::Type{DynmatNamelist})
    asr = ASR[request(
        io,
        @green("Please select the type of acoustic sum rule `asr`: "),
        RadioMenu(ASR; charset=:ascii),
    )]
    print(io, @green "Please input mass for each atom type `amass` (separated by spaces): ")
    amass = map(Base.Fix1(parse, Float64), split(readline(io), " "; keepempty=false))
    dynmat = DynmatNamelist(; asr, amass)
    return FieldSetter()(io, dynmat)
end
(builder::InputBuilder)(io::IO, ::Type{PhInput}) = PhInput(builder(io, PhNamelist))
(builder::InputBuilder)(io::IO, ::Type{Q2rInput}) = Q2rInput(builder(io, Q2rNamelist))
(builder::InputBuilder)(io::IO, ::Type{MatdynInput}) = MatdynInput(builder(io, MatdynNamelist))
(builder::InputBuilder)(io::IO, ::Type{DynmatInput}) = DynmatInput(builder(io, DynmatNamelist))

end
