module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: RED_FG, GREEN_FG
using QuantumESPRESSOBase: entryname, qestring
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

import ..QuantumESPRESSOHelpers: build, setfield_helper

function build(terminal::TTYTerminal, ::Type{ControlNamelist})
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    restart_modes = pairs(("from_scratch", "restart"))
    calculation = calculations[request(
        terminal,
        string(GREEN_FG("What exact calculation do you want to run?")),
        RadioMenu(collect(values(calculations))),
    )]
    restart_mode = restart_modes[request(
        terminal, string(GREEN_FG("Starting from scratch?")), RadioMenu(["yes", "no"])
    )]
    print(terminal, string(GREEN_FG("Convergence threshold on total energy (a.u): ")))
    etot_conv_thr = parse(Float64, readline(terminal))
    print(terminal, string(GREEN_FG("Convergence threshold on forces (a.u): ")))
    forc_conv_thr = parse(Float64, readline(terminal))
    control = T(;
        calculation=calculation,
        restart_mode=restart_mode,
        etot_conv_thr=etot_conv_thr,
        forc_conv_thr=forc_conv_thr,
    )
    return setfield_helper(terminal, control)
end
function build(terminal::TTYTerminal, ::Type{SystemNamelist})
    print(terminal, string(GREEN_FG("Please input the Bravais lattice index `ibrav`: ")))
    ibrav = parse(Int, readline(terminal))
    print(terminal, string(GREEN_FG("Please input a `celldm` 1-6 (separated by spaces): ")))
    celldm = map(x -> parse(Float64, x), split(readline(terminal), " "; keepempty=false))
    print(
        terminal,
        string(GREEN_FG("Please input the number of atoms in the unit cell `nat`: ")),
    )
    nat = parse(Int, readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG("Please input the number of types of atoms in the unit cell `ntyp`: ")
        ),
    )
    ntyp = parse(Int, readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: "
            ),
        ),
    )
    ecutwfc = parse(Float64, readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: "
            ),
        ),
    )
    ecutrho = parse(Float64, readline(terminal))
    system = T(;
        ibrav=ibrav, celldm=celldm, nat=nat, ntyp=ntyp, ecutwfc=ecutwfc, ecutrho=ecutrho
    )
    return setfield_helper(terminal, system)
end
function build(terminal::TTYTerminal, ::Type{ElectronsNamelist})
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the convergence threshold for selfconsistency `conv_thr`: "
            ),
        ),
    )
    conv_thr = parse(Float64, readline(terminal))
    diagonalizations = pairs(("david", "cg", "cg-serial", "david-serial"))
    diagonalization = diagonalizations[request(
        terminal,
        string(GREEN_FG("Please input the diagonalization method `diagonalization`: ")),
        RadioMenu(collect(values(diagonalizations))),
    )]
    electrons = T(; conv_thr=conv_thr, diagonalization=diagonalization)
    return setfield_helper(terminal, electrons)
end
function build(terminal::TTYTerminal, ::Type{IonsNamelist})
    ion_dynamics_pool = pairs((
        "none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"
    ))
    ion_dynamics = ion_dynamics_pool[request(
        terminal,
        string(GREEN_FG("Please input the type of ionic dynamics `ion_dynamics`: ")),
        RadioMenu(collect(values(ion_dynamics_pool))),
    )]
    ion_temperature_pool = (
        "rescaling",
        "rescale-v",
        "rescale-T",
        "reduce-T",
        "berendsen",
        "andersen",
        "initial",
        "not_controlled",
    )
    ion_temperature = ion_temperature_pool[request(
        terminal,
        string(GREEN_FG("Please input the ions temperature `ion_temperature`: ")),
        RadioMenu(collect(values(ion_temperature_pool))),
    )]
    ions = T(; ion_dynamics=ion_dynamics, ion_temperature=ion_temperature)
    return setfield_helper(terminal, ions)
end
function build(terminal::TTYTerminal, ::Type{CellNamelist})
    cell_dynamics_pool = pairs(("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"))
    cell_dynamics = cell_dynamics_pool[request(
        terminal,
        string(
            GREEN_FG("Please input the type of dynamics for the cell `cell_dynamics`: ")
        ),
        RadioMenu(collect(values(cell_dynamics_pool))),
    )]
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: ",
            ),
        ),
    )
    press = parse(Float64, readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: ",
            ),
        ),
    )
    wmass = parse(Float64, readline(terminal))
    print(
        terminal,
        string(
            GREEN_FG(
                "Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: ",
            ),
        ),
    )
    press_conv_thr = parse(Float64, readline(terminal))
    cell = T(;
        cell_dynamics=cell_dynamics, press=press, wmass=wmass, press_conv_thr=press_conv_thr
    )
    return setfield_helper(terminal, cell)
end

function build(terminal::TTYTerminal, ::Type{T}) where {T<:KPointsCard}
    kpt_style = request(
        terminal,
        string(GREEN_FG("What k-point style do you want?")),
        RadioMenu(["gamma", "automatic"]),
    )
    return if kpt_style == 1
        KPointsCard("gamma", GammaPoint())
    else  # "automatic"
        print(
            terminal,
            string(
                GREEN_FG("What 3-element k-point grid do you want (separated by spaces): ")
            ),
        )
        grid = map(x -> parse(Int, x), split(readline(terminal), " "; keepempty=false))
        print(
            terminal,
            string(
                GREEN_FG(
                    "What 3-element k-point offsets do you want (separated by spaces): "
                ),
            ),
        )
        offsets = map(x -> parse(Int, x), split(readline(terminal), " "; keepempty=false))
        return KPointsCard("automatic", MonkhorstPackGrid(grid, offsets))
    end
end

function build(terminal::TTYTerminal, ::Type{PWInput})
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, entryname(S) => build(terminal, S))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, string(RED_FG("Something wrong happens, try again!")))
            end
        end
    end
    if fields[:control].calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, entryname(IonsNamelist) => build(terminal, IonsNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, string(RED_FG("Something wrong happens, try again!")))
            end
        end
    else
        IonsNamelist()
    end
    if fields[:control].calculation ∈ ("vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(fields, entryname(CellNamelist) => build(terminal, CellNamelist))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, string(RED_FG("Something wrong happens, try again!")))
            end
        end
    else
        CellNamelist()
    end
    haserror = true
    while haserror
        try
            push!(fields, entryname(KPointsCard) => build(terminal, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(terminal, string(RED_FG("Something wrong happens, try again!")))
        end
    end
    push!(fields, entryname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(
        fields,
        entryname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]),
    )
    result = T(; fields...)
    saveresult = pairs((true, false))[request(
        terminal,
        string(GREEN_FG("Do you want to save the generated input to file?")),
        RadioMenu(["yes", "no"]),
    )]
    if saveresult
        print(terminal, string(GREEN_FG("Input file name: ")))
        write(chomp(readline(terminal)), qestring(result))
    end
    return result
end

end
