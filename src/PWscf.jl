module PWscf

using REPL.TerminalMenus: RadioMenu, request

using AbInitioSoftwareBase: groupname
using QuantumESPRESSOBase.PWscf:
    ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    AtomicSpecies,
    AtomicSpeciesCard,
    AtomicPosition,
    AtomicPositionsCard,
    KPointsCard,
    CellParametersCard,
    GammaPoint,
    MonkhorstPackGrid,
    KPointsCard,
    PWInput
using REPL.TerminalMenus: RadioMenu, request, terminal
using Term: @green, @red

import ..QuantumESPRESSOHelpers: build, setfield_helper

const CALCULATIONS = Base.vect("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
const RESTART_MODES = Base.vect("from_scratch", "restart")
const DIAGONALIZATIONS = Base.vect("david", "cg", "cg-serial", "david-serial")
const ION_DYNAMICS_POOL = Base.vect(
    "none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"
)
const ION_TEMPERATURES = Base.vect(
    "rescaling",
    "rescale-v",
    "rescale-T",
    "reduce-T",
    "berendsen",
    "andersen",
    "initial",
    "not_controlled",
)
const CELL_DYNAMICS_POOL = Base.vect("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w")

function build(term::IO, ::Type{ControlNamelist})
    calculation = CALCULATIONS[request(
        term,
        @green("What exact calculation do you want to run?"),
        RadioMenu(CALCULATIONS; charset=:ascii),
    )]
    restart_mode = RESTART_MODES[request(
        term, @green("Starting from scratch?"), RadioMenu(["yes", "no"])
    )]
    print(term, @green "Convergence threshold on total energy (a.u): ")
    etot_conv_thr = parse(Float64, readline(term))
    print(term, @green "Convergence threshold on forces (a.u): ")
    forc_conv_thr = parse(Float64, readline(term))
    control = ControlNamelist(; calculation, restart_mode, etot_conv_thr, forc_conv_thr)
    return setfield_helper(term, control)
end
function build(term::IO, ::Type{SystemNamelist})
    print(term, @green "Please input the Bravais lattice index `ibrav`: ")
    ibrav = parse(Int, readline(term))
    print(term, @green "Please input a `celldm` 1-6 (separated by spaces): ")
    celldm = map(Base.Fix1(parse, Float64), split(readline(term), " "; keepempty=false))
    print(term, @green "Please input the number of atoms in the unit cell `nat`: ")
    nat = parse(Int, readline(term))
    print(
        term, @green "Please input the number of types of atoms in the unit cell `ntyp`: "
    )
    ntyp = parse(Int, readline(term))
    print(
        term,
        @green "Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: "
    )
    ecutwfc = parse(Float64, readline(term))
    print(
        term,
        @green "Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: "
    )
    ecutrho = parse(Float64, readline(term))
    system = SystemNamelist(; ibrav, celldm, nat, ntyp, ecutwfc, ecutrho)
    return setfield_helper(term, system)
end
function build(term::IO, ::Type{ElectronsNamelist})
    print(
        term,
        @green "Please input the convergence threshold for selfconsistency `conv_thr`: "
    )
    conv_thr = parse(Float64, readline(term))
    diagonalization = DIAGONALIZATIONS[request(
        term,
        @green("Please input the diagonalization method `diagonalization`: "),
        RadioMenu(DIAGONALIZATIONS; charset=:ascii),
    )]
    electrons = ElectronsNamelist(; conv_thr, diagonalization)
    return setfield_helper(term, electrons)
end
function build(term::IO, ::Type{IonsNamelist})
    ion_dynamics = ION_DYNAMICS_POOL[request(
        term,
        @green("Please input the type of ionic dynamics `ion_dynamics`: "),
        RadioMenu(ION_DYNAMICS_POOL; charset=:ascii),
    )]
    ion_temperature = ION_TEMPERATURES[request(
        term,
        @green("Please input the ions temperature `ion_temperature`: "),
        RadioMenu(ION_TEMPERATURES; charset=:ascii),
    )]
    ions = IonsNamelist(; ion_dynamics, ion_temperature)
    return setfield_helper(term, ions)
end
function build(term::IO, ::Type{CellNamelist})
    cell_dynamics = CELL_DYNAMICS_POOL[request(
        term,
        @green("Please input the type of dynamics for the cell `cell_dynamics`: "),
        RadioMenu(CELL_DYNAMICS_POOL; charset=:ascii),
    )]
    print(
        term,
        @green "Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: "
    )
    press = parse(Float64, readline(term))
    print(
        term,
        @green "Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: "
    )
    wmass = parse(Float64, readline(term))
    print(
        term,
        @green "Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: "
    )
    press_conv_thr = parse(Float64, readline(term))
    cell = CellNamelist(; cell_dynamics, press, wmass, press_conv_thr)
    return setfield_helper(term, cell)
end
function build(term::IO, ::Type{KPointsCard})
    kpt_style = request(
        term, @green("What k-point style do you want?"), RadioMenu(["gamma", "automatic"])
    )
    return if kpt_style == 1
        GammaPointCard()
    else  # "automatic"
        print(
            term, @green "What 3-element k-point grid do you want (separated by spaces): "
        )
        grid = map(x -> parse(Int, x), split(readline(term), " "; keepempty=false))
        print(
            term,
            @green "What 3-element k-point offsets do you want (separated by spaces): "
        )
        offsets = map(x -> parse(Int, x), split(readline(term), " "; keepempty=false))
        return KMeshCard(MonkhorstPackGrid(grid, offsets))
    end
end
function build(term::IO, ::Type{PWInput})
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, groupname(S) => build(term, S))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(term, @red "Something wrong happens, try again!")
            end
        end
    end
    if fields[:control].calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, groupname(IonsNamelist) => build(term, IonsNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(term, @red "Something wrong happens, try again!")
            end
        end
    else
        IonsNamelist()
    end
    if fields[:control].calculation ∈ ("vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, groupname(CellNamelist) => build(term, CellNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(term, @red "Something wrong happens, try again!")
            end
        end
    else
        CellNamelist()
    end
    haserror = true
    while haserror
        try
            push!(fields, groupname(KPointsCard) => build(term, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(term, @red "Something wrong happens, try again!")
        end
    end
    push!(fields, groupname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(
        fields,
        groupname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]),
    )
    result = PWInput(; fields...)
    saveresult = Base.vect(true, false)[request(
        term,
        @green("Do you want to save the generated input to file?"),
        RadioMenu(["yes", "no"]),
    )]
    if saveresult
        print(term, @green "Input file name: ")
        write(chomp(readline(term)), result)
    end
    return result
end
build(T::Type) = build(terminal, T)

end
