"""
RheoJL

A Julia package for rheology-related computations.
"""
module RheoJL

    using CSV, DataFrames, Interpolations

    include("utils.jl")
    include("squeeze_flow.jl")
    include("solver.jl")
    include("rheology.jl")

    # Exports from "utils.jl"
    export vertex_to_midpoint, right_integrate, integral, list_data_files, load_data, geometric_time_sequence

    # Exports from "squeeze_flow.jl"
    export shear_stress, velocity_profile, flux, force

    # Exports from "solver.jl"
    export solve_∂p∂r, solve_system, integrate_system, newton, NewtonDidNotConverge

    # Exports from "rheology.jl"
    export solve_τ

end # module RheoJL