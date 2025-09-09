using CSV, DataFrames

"""
    right_integrate(rᵥ, fₘ)

Compute the right integral of a function sampled at midpoints.

Given vectors `rᵥ` (vertices, length N) and `fₘ` (midpoints, length N-1), returns a vector of length N where each entry is the right integral from the end, assuming zero at the rightmost vertex.

# Examples
```jldoctest
julia> using RheoJL

julia> rᵥ = [0.0, 1.0, 3.0]; fₘ = [2.0, 4.0];

julia> right_integrate(rᵥ, fₘ)
3-element Vector{Float64}:
 -10.0
  -8.0
   0.0
```
"""
function right_integrate(rᵥ, fₘ)

    Fᵥ = similar(rᵥ)
    Fᵥ[end] = 0.0

    for i in length(fₘ):-1:1
        Δr = rᵥ[i+1] - rᵥ[i]
        Fᵥ[i] = Fᵥ[i+1] - fₘ[i]*Δr
    end

    return Fᵥ
end

"""
    integral(rᵥ, fₘ)

Compute the definite integral of a function sampled at midpoints.

Given vectors `rᵥ` (vertices, length N) and `fₘ` (midpoints, length N-1), returns the sum over all intervals.

# Examples
```jldoctest
julia> using RheoJL

julia> rᵥ = [0.0, 1.0, 3.0]; fₘ = [2.0, 4.0];

julia> integral(rᵥ, fₘ)
10.0
```
"""
function integral(rᵥ, fₘ)
    F = 0.0
    for i in 1:length(rᵥ)-1
        Δr = rᵥ[i+1] - rᵥ[i]
        F += fₘ[i]*Δr
    end

    return F
end

"""
    vertex_to_midpoint(fᵥ)

Convert a vector of values defined at vertices to values at midpoints between vertices.

Given a vector `fᵥ` of length N, returns a vector of length N-1 where each element is the average of two consecutive elements of `fᵥ`.

# Examples
```jldoctest
julia> using RheoJL

julia> vertex_to_midpoint([1.0, 3.0, 7.0])
2-element Vector{Float64}:
 2.0
 5.0
```
"""
function vertex_to_midpoint(fᵥ)
    return 0.5*(fᵥ[1:end-1]+fᵥ[2:end])
end

function load_data(name)
    data_path = joinpath(@__DIR__, "..", "data", name)
    return CSV.read(data_path, DataFrame)
end