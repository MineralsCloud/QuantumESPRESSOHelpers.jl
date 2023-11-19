using AbInitioSoftwareBase: groupname
using CrystallographyBase: MonkhorstPackGrid
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
    GammaPointCard,
    KMeshCard,
    PWInput
using QuantumESPRESSOFormatter.PWscf
using REPL.TerminalMenus: RadioMenu, request
using Term: @green, @red

using ..QuantumESPRESSOHelpers: help_set

const CALCULATION = Base.vect("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md")
const RESTART_MODE = Base.vect("from_scratch", "restart")
const DIAGONALIZATION = Base.vect("david", "cg", "cg-serial", "david-serial")
const ION_DYNAMICS = Base.vect(
    "none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"
)
const ION_TEMPERATURE = Base.vect(
    "rescaling",
    "rescale-v",
    "rescale-T",
    "reduce-T",
    "berendsen",
    "andersen",
    "initial",
    "not_controlled",
)
const CELL_DYNAMICS = Base.vect("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w")

function build(io::IO, ::Type{ControlNamelist})
    calculation = CALCULATION[request(
        io,
        @green("What exact calculation do you want to run?"),
        RadioMenu(CALCULATION; charset=:ascii),
    )]
    restart_mode = RESTART_MODE[request(
        io, @green("Starting from scratch?"), RadioMenu(["yes", "no"])
    )]
    print(io, @green "Convergence threshold on total energy (a.u): ")
    etot_conv_thr = parse(Float64, readline(io))
    print(io, @green "Convergence threshold on forces (a.u): ")
    forc_conv_thr = parse(Float64, readline(io))
    control = ControlNamelist(; calculation, restart_mode, etot_conv_thr, forc_conv_thr)
    return help_set(io, control)
end
function build(io::IO, ::Type{SystemNamelist})
    print(io, @green "Please input the Bravais lattice index `ibrav`: ")
    ibrav = parse(Int, readline(io))
    print(io, @green "Please input a `celldm` 1-6 (separated by spaces): ")
    celldm = map(Base.Fix1(parse, Float64), split(readline(io), " "; keepempty=false))
    print(io, @green "Please input the number of atoms in the unit cell `nat`: ")
    nat = parse(Int, readline(io))
    print(io, @green "Please input the number of types of atoms in the unit cell `ntyp`: ")
    ntyp = parse(Int, readline(io))
    print(
        io,
        @green "Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: "
    )
    ecutwfc = parse(Float64, readline(io))
    print(
        io,
        @green "Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: "
    )
    ecutrho = parse(Float64, readline(io))
    system = SystemNamelist(; ibrav, celldm, nat, ntyp, ecutwfc, ecutrho)
    return help_set(io, system)
end
function build(io::IO, ::Type{ElectronsNamelist})
    print(
        io, @green "Please input the convergence threshold for selfconsistency `conv_thr`: "
    )
    conv_thr = parse(Float64, readline(io))
    diagonalization = DIAGONALIZATION[request(
        io,
        @green("Please input the diagonalization method `diagonalization`: "),
        RadioMenu(DIAGONALIZATION; charset=:ascii),
    )]
    electrons = ElectronsNamelist(; conv_thr, diagonalization)
    return help_set(io, electrons)
end
function build(io::IO, ::Type{IonsNamelist})
    ion_dynamics = ION_DYNAMICS[request(
        io,
        @green("Please input the type of ionic dynamics `ion_dynamics`: "),
        RadioMenu(ION_DYNAMICS; charset=:ascii),
    )]
    ion_temperature = ION_TEMPERATURE[request(
        io,
        @green("Please input the ions temperature `ion_temperature`: "),
        RadioMenu(ION_TEMPERATURE; charset=:ascii),
    )]
    ions = IonsNamelist(; ion_dynamics, ion_temperature)
    return help_set(io, ions)
end
function build(io::IO, ::Type{CellNamelist})
    cell_dynamics = CELL_DYNAMICS[request(
        io,
        @green("Please input the type of dynamics for the cell `cell_dynamics`: "),
        RadioMenu(CELL_DYNAMICS; charset=:ascii),
    )]
    print(
        io,
        @green "Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: "
    )
    press = parse(Float64, readline(io))
    print(
        io,
        @green "Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: "
    )
    wmass = parse(Float64, readline(io))
    print(
        io,
        @green "Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: "
    )
    press_conv_thr = parse(Float64, readline(io))
    cell = CellNamelist(; cell_dynamics, press, wmass, press_conv_thr)
    return help_set(io, cell)
end
function build(io::IO, ::Type{KPointsCard})
    kpt_style = request(
        io, @green("What k-point style do you want?"), RadioMenu(["gamma", "automatic"])
    )
    return if kpt_style == 1
        GammaPointCard()
    else  # "automatic"
        print(io, @green "What 3-element k-point grid do you want (separated by spaces): ")
        grid = map(Base.Fix1(parse, Int), split(readline(io), " "; keepempty=false))
        print(
            io, @green "What 3-element k-point offsets do you want (separated by spaces): "
        )
        offsets = map(Base.Fix1(parse, Int), split(readline(io), " "; keepempty=false))
        return KMeshCard(MonkhorstPackGrid(grid, offsets))
    end
end
function build(io::IO, ::Type{PWInput})
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, groupname(S) => build(io, S))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(io, @red "Something wrong happens, try again!")
            end
        end
    end
    if fields[:control].calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, groupname(IonsNamelist) => build(io, IonsNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(io, @red "Something wrong happens, try again!")
            end
        end
    else
        IonsNamelist()
    end
    if fields[:control].calculation ∈ ("vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, groupname(CellNamelist) => build(io, CellNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(io, @red "Something wrong happens, try again!")
            end
        end
    else
        CellNamelist()
    end
    haserror = true
    while haserror
        try
            push!(fields, groupname(KPointsCard) => build(io, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(io, @red "Something wrong happens, try again!")
        end
    end
    push!(fields, groupname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(
        fields,
        groupname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]),
    )
    result = PWInput(; fields...)
    saveresult = Base.vect(true, false)[request(
        io, @green("Do you want to save the generated input to file?"), YES_NO_MENU
    )]
    if saveresult
        print(io, @green "Input file name: ")
        write(chomp(readline(io)), result)
    end
    return result
end
