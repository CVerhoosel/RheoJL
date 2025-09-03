"""
RheoJL

A Julia package for rheology-related computations.
"""
module RheoJL

    include("utils.jl")
    include("squeeze_flow.jl")
    include("solver.jl")

    # Exports from "utils.jl"
    export vertex_to_midpoint, right_integrate, integral
    
    # Exports from "squeeze_flow.jl"
    export velocity_profile, flux

    # Exports from "solver.jl"
    export solve_∂p∂r

end # module RheoJL