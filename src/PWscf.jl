module PWscf

using REPL.TerminalMenus: RadioMenu, request

using Crayons.Box: RED_FG, GREEN_FG
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

import ..QuantumESPRESSOHelpers: build, setfield_helper

function build(term::IO, ::Type{ControlNamelist})
    calculations = pairs(("scf", "nscf", "bands", "relax", "md", "vc-relax", "vc-md"))
    restart_modes = pairs(("from_scratch", "restart"))
    calculation = calculations[request(
        term,
        string(GREEN_FG("What exact calculation do you want to run?")),
        RadioMenu(collect(values(calculations))),
    )]
    restart_mode = restart_modes[request(
        term, string(GREEN_FG("Starting from scratch?")), RadioMenu(["yes", "no"])
    )]
    print(term, string(GREEN_FG("Convergence threshold on total energy (a.u): ")))
    etot_conv_thr = parse(Float64, readline(term))
    print(term, string(GREEN_FG("Convergence threshold on forces (a.u): ")))
    forc_conv_thr = parse(Float64, readline(term))
    control = T(;
        calculation=calculation,
        restart_mode=restart_mode,
        etot_conv_thr=etot_conv_thr,
        forc_conv_thr=forc_conv_thr,
    )
    return setfield_helper(term, control)
end
function build(term::IO, ::Type{SystemNamelist})
    print(term, string(GREEN_FG("Please input the Bravais lattice index `ibrav`: ")))
    ibrav = parse(Int, readline(term))
    print(term, string(GREEN_FG("Please input a `celldm` 1-6 (separated by spaces): ")))
    celldm = map(x -> parse(Float64, x), split(readline(term), " "; keepempty=false))
    print(
        term, string(GREEN_FG("Please input the number of atoms in the unit cell `nat`: "))
    )
    nat = parse(Int, readline(term))
    print(
        term,
        string(
            GREEN_FG("Please input the number of types of atoms in the unit cell `ntyp`: ")
        ),
    )
    ntyp = parse(Int, readline(term))
    print(
        term,
        string(
            GREEN_FG(
                "Please input the kinetic energy cutoff (Ry) for wavefunctions `ecutwfc`: "
            ),
        ),
    )
    ecutwfc = parse(Float64, readline(term))
    print(
        term,
        string(
            GREEN_FG(
                "Please input the Kinetic energy cutoff (Ry) for charge density `ecutrho`: "
            ),
        ),
    )
    ecutrho = parse(Float64, readline(term))
    system = T(;
        ibrav=ibrav, celldm=celldm, nat=nat, ntyp=ntyp, ecutwfc=ecutwfc, ecutrho=ecutrho
    )
    return setfield_helper(term, system)
end
function build(term::IO, ::Type{ElectronsNamelist})
    print(
        term,
        string(
            GREEN_FG(
                "Please input the convergence threshold for selfconsistency `conv_thr`: "
            ),
        ),
    )
    conv_thr = parse(Float64, readline(term))
    diagonalizations = pairs(("david", "cg", "cg-serial", "david-serial"))
    diagonalization = diagonalizations[request(
        term,
        string(GREEN_FG("Please input the diagonalization method `diagonalization`: ")),
        RadioMenu(collect(values(diagonalizations))),
    )]
    electrons = T(; conv_thr=conv_thr, diagonalization=diagonalization)
    return setfield_helper(term, electrons)
end
function build(term::IO, ::Type{IonsNamelist})
    ion_dynamics_pool = pairs((
        "none", "bfgs", "damp", "verlet", "langevin", "langevin-smc", "beeman"
    ))
    ion_dynamics = ion_dynamics_pool[request(
        term,
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
        term,
        string(GREEN_FG("Please input the ions temperature `ion_temperature`: ")),
        RadioMenu(collect(values(ion_temperature_pool))),
    )]
    ions = T(; ion_dynamics=ion_dynamics, ion_temperature=ion_temperature)
    return setfield_helper(term, ions)
end
function build(term::IO, ::Type{CellNamelist})
    cell_dynamics_pool = pairs(("none", "sd", "damp-pr", "damp-w", "bfgs", "pr", "w"))
    cell_dynamics = cell_dynamics_pool[request(
        term,
        string(
            GREEN_FG("Please input the type of dynamics for the cell `cell_dynamics`: ")
        ),
        RadioMenu(collect(values(cell_dynamics_pool))),
    )]
    print(
        term,
        string(
            GREEN_FG(
                "Please input the target pressure [KBar] in a variable-cell md or relaxation run `press`: ",
            ),
        ),
    )
    press = parse(Float64, readline(term))
    print(
        term,
        string(
            GREEN_FG(
                "Please input the fictitious cell mass [amu] for variable-cell simulations `wmass`: ",
            ),
        ),
    )
    wmass = parse(Float64, readline(term))
    print(
        term,
        string(
            GREEN_FG(
                "Please input the Convergence threshold on the pressure for variable cell `press_conv_thr`: ",
            ),
        ),
    )
    press_conv_thr = parse(Float64, readline(term))
    cell = T(;
        cell_dynamics=cell_dynamics, press=press, wmass=wmass, press_conv_thr=press_conv_thr
    )
    return setfield_helper(term, cell)
end

function build(term::IO, ::Type{T}) where {T<:KPointsCard}
    kpt_style = request(
        term,
        string(GREEN_FG("What k-point style do you want?")),
        RadioMenu(["gamma", "automatic"]),
    )
    return if kpt_style == 1
        KPointsCard("gamma", GammaPoint())
    else  # "automatic"
        print(
            term,
            string(
                GREEN_FG("What 3-element k-point grid do you want (separated by spaces): ")
            ),
        )
        grid = map(x -> parse(Int, x), split(readline(term), " "; keepempty=false))
        print(
            term,
            string(
                GREEN_FG(
                    "What 3-element k-point offsets do you want (separated by spaces): "
                ),
            ),
        )
        offsets = map(x -> parse(Int, x), split(readline(term), " "; keepempty=false))
        return KPointsCard("automatic", MonkhorstPackGrid(grid, offsets))
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
                println(term, string(RED_FG("Something wrong happens, try again!")))
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
                println(term, string(RED_FG("Something wrong happens, try again!")))
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
                println(term, string(RED_FG("Something wrong happens, try again!")))
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
            println(term, string(RED_FG("Something wrong happens, try again!")))
        end
    end
    push!(fields, groupname(AtomicSpeciesCard) => AtomicSpeciesCard(AtomicSpecies[]))
    push!(
        fields,
        groupname(AtomicPositionsCard) => AtomicPositionsCard("alat", AtomicPosition[]),
    )
    result = T(; fields...)
    saveresult = pairs((true, false))[request(
        term,
        string(GREEN_FG("Do you want to save the generated input to file?")),
        RadioMenu(["yes", "no"]),
    )]
    if saveresult
        print(term, string(GREEN_FG("Input file name: ")))
        write(chomp(readline(term)), result)
    end
    return result
end

end
