module PWscf

using REPL.Terminals: TTYTerminal
using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: GREEN_FG
using QuantumESPRESSOBase.Inputs.PWscf:
    ControlNamelist, SystemNamelist, ElectronsNamelist, IonsNamelist, CellNamelist

using ..Namelists

function Namelists.namelist_builder(
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
    return Namelists.setfield_helper(terminal, control)
end # function namelist_builder
function Namelists.namelist_builder(
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
    return Namelists.setfield_helper(terminal, system)
end # function namelist_builder
function Namelists.namelist_builder(
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
    return Namelists.setfield_helper(terminal, electrons)
end # function namelist_builder
function Namelists.namelist_builder(
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
    return Namelists.setfield_helper(terminal, ions)
end # function namelist_builder
function Namelists.namelist_builder(
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
    return Namelists.setfield_helper(terminal, cell)
end # function namelist_builder

end # module PWscf
