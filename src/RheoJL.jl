"""
RheoJL

A Julia package for rheology-related computations.
"""
module RheoJL

    using ProgressMeter

    include("utils.jl")
    include("squeeze_flow.jl")
    include("solver.jl")

    # Exports from "utils.jl"
    export vertex_to_midpoint, right_integrate, integral, load_data

    # Exports from "squeeze_flow.jl"
    export velocity_profile, flux, force

    # Exports from "solver.jl"
    export solve_∂p∂r, solve_system, integrate_system,newton

end # module RheoJL