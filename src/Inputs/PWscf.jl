module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: RED_FG, GREEN_FG
using QuantumESPRESSOBase.Inputs: entryname, qestring
using QuantumESPRESSOBase.Inputs.PWscf:
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
    GammaPoint, MonkhorstPackGrid, KPointsCard,
    PWInput

using ..Inputs: namelist_builder, card_builder, input_builder, setfield_helper

import QuantumESPRESSOHelpers.Inputs

function Inputs.namelist_builder(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:ControlNamelist}
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    restart_modes = pairs(("from_scratch", "restart"))
    calculation = calculations[request(
        terminal,
        GREEN_FG("What exact calculation do you want to run?") |> string,
        RadioMenu(collect(values(calculations))),
    )]
    restart_mode = restart_modes[request(
        terminal,
        GREEN_FG("Starting from scratch?") |> string,
        RadioMenu(["yes", "no"]),
    )]
    print(terminal, GREEN_FG("Convergence threshold on total energy (a.u): ") |> string)
    etot_conv_thr = parse(Float64, readline(terminal))
    print(terminal, GREEN_FG("Convergence threshold on forces (a.u): ") |> string)
    forc_conv_thr = parse(Float64, readline(terminal))
    control = T(
        calculation = calculation,
        restart_mode = restart_mode,
        etot_conv_thr = etot_conv_thr,
        forc_conv_thr = forc_conv_thr,
    )
    return setfield_helper(terminal, control)
end # function namelist_builder
function Inputs.namelist_builder(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:SystemNamelist}
    print(terminal, GREEN_FG("Please input the Bravais lattice index `ibrav`: ") |> string)
    ibrav = parse(Int, readline(terminal))
    print(
        terminal,
        GREEN_FG("Please input a `celldm` 1-6 (separated by spaces): ") |> string,
    )
    celldm = map(x -> parse(Float64, x), split(readline(terminal), " ", keepempty = false))
    print(
        terminal,
        GREEN_FG("Please input the number of atoms in the unit cell `nat`: ") |> string,
    )
    nat = parse(Int, readline(terminal))
    print(
        terminal,
        GREEN_FG("Please input the number of types of atoms in the unit cell `ntyp`: ") |>
        string,
    )
    ntyp = parse(Int, readline(terminal))
    print(
        terminal,
        GREEN_FG("Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: ") |>
        string,
    )
    ecutwfc = parse(Float64, readline(terminal))
    print(
        terminal,
        GREEN_FG("Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: ") |>
        string,
    )
    ecutrho = parse(Float64, readline(terminal))
    system = T(
        ibrav = ibrav,
        celldm = celldm,
        nat = nat,
        ntyp = ntyp,
        ecutwfc = ecutwfc,
        ecutrho = ecutrho,
    )
    return setfield_helper(terminal, system)
end # function namelist_builder
function Inputs.namelist_builder(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:ElectronsNamelist}
    print(
        terminal,
        GREEN_FG("Please input the convergence threshold for selfconsistency `conv_thr`: ") |>
        string,
    )
    conv_thr = parse(Float64, readline(terminal))
    diagonalizations = pairs(("david", "cg", "cg-serial", "david-serial"))
    diagonalization = diagonalizations[request(
        terminal,
        GREEN_FG("Please input the diagonalization method `diagonalization`: ") |> string,
        RadioMenu(collect(values(diagonalizations))),
    )]
    electrons = T(conv_thr = conv_thr, diagonalization = diagonalization)
    return setfield_helper(terminal, electrons)
end # function namelist_builder
function Inputs.namelist_builder(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:IonsNamelist}
    ion_dynamics_pool =
        pairs(("none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"))
    ion_dynamics = ion_dynamics_pool[request(
        terminal,
        GREEN_FG("Please input the type of ionic dynamics `ion_dynamics`: ") |> string,
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
        GREEN_FG("Please input the ions temperature `ion_temperature`: ") |> string,
        RadioMenu(collect(values(ion_temperature_pool))),
    )]
    ions = T(ion_dynamics = ion_dynamics, ion_temperature = ion_temperature)
    return setfield_helper(terminal, ions)
end # function namelist_builder
function Inputs.namelist_builder(
    terminal::TTYTerminal,
    ::Type{T},
) where {T<:CellNamelist}
    cell_dynamics_pool = pairs(("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"))
    cell_dynamics = cell_dynamics_pool[request(
        terminal,
        GREEN_FG("Please input the type of dynamics for the cell `cell_dynamics`: ") |>
        string,
        RadioMenu(collect(values(cell_dynamics_pool))),
    )]
    print(
        terminal,
        GREEN_FG("Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: ") |>
        string,
    )
    press = parse(Float64, readline(terminal))
    print(
        terminal,
        GREEN_FG("Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: ") |>
        string,
    )
    wmass = parse(Float64, readline(terminal))
    print(
        terminal,
        GREEN_FG("Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: ") |>
        string,
    )
    press_conv_thr = parse(Float64, readline(terminal))
    cell = T(
        cell_dynamics = cell_dynamics,
        press = press,
        wmass = wmass,
        press_conv_thr = press_conv_thr,
    )
    return setfield_helper(terminal, cell)
end # function namelist_builder

function Inputs.card_builder(terminal::TTYTerminal, ::Type{T}) where {T<:KPointsCard}
    kpt_style = request(
        terminal,
        GREEN_FG("What k-point style do you want?") |> string,
        RadioMenu(["gamma", "automatic"]),
    )
    return if kpt_style == 1
        KPointsCard("gamma", GammaPoint())
    else  # "automatic"
        print(
            terminal,
            GREEN_FG("What 3-element k-point grid do you want (separated by spaces): ") |>
            string,
        )
        grid = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        print(
            terminal,
            GREEN_FG("What 3-element k-point offsets do you want (separated by spaces): ") |>
            string,
        )
        offsets = map(x -> parse(Int, x), split(readline(terminal), " ", keepempty = false))
        return KPointsCard("automatic", MonkhorstPackGrid(grid, offsets))
    end
end # function card_builder

function Inputs.input_builder(terminal::TTYTerminal, ::Type{T}) where {T<:PWInput}
    fields = Dict{Symbol,Any}()
    for S in (ControlNamelist, SystemNamelist, ElectronsNamelist)
        haserror = true
        while haserror
            try
                push!(fields, entryname(S) => namelist_builder(terminal, S))
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, RED_FG("Something wrong happens, try again!") |> string)
            end
        end
    end
    if fields[:control].calculation ∈ ("relax", "md", "vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(
                    fields,
                    entryname(IonsNamelist) => namelist_builder(terminal, IonsNamelist),
                )
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, RED_FG("Something wrong happens, try again!") |> string)
            end
        end
    else
        IonsNamelist()
    end
    if fields[:control].calculation ∈ ("vc-relax", "vc-md")
        haserror = true
        while haserror
            try
                push!(
                    fields,
                    entryname(CellNamelist) => namelist_builder(terminal, CellNamelist),
                )
                haserror = false
            catch e
                isa(e, InterruptException) && rethrow(e)
                println(terminal, RED_FG("Something wrong happens, try again!") |> string)
            end
        end
    else
        CellNamelist()
    end
    haserror = true
    while haserror
        try
            push!(fields, entryname(KPointsCard) => card_builder(terminal, KPointsCard))
            haserror = false
        catch e
            isa(e, InterruptException) && rethrow(e)
            println(terminal, RED_FG("Something wrong happens, try again!") |> string)
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
        GREEN_FG("Do you want to save the generated input to file?") |> string,
        RadioMenu(["yes", "no"]),
    )]
    if saveresult
        print(terminal, GREEN_FG("Input file name: ") |> string)
        write(readline(terminal) |> chomp, qestring(result))
    end
    return result
end # function input_builder

end # module PWscf
