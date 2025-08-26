function right_integrate(rᵥ, fₘ)

    Fᵥ = similar(rᵥ)
    Fᵥ[end] = 0.0

    for i in length(fₘ):-1:1
        Δr = rᵥ[i+1] - rᵥ[i]
        Fᵥ[i] = Fᵥ[i+1] - fₘ[i]*Δr
    end

    return Fᵥ
end

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
julia> vertex_to_midpoint([1.0, 3.0, 7.0])
2-element Vector{Float64}:
 2.0
 5.0
```
"""
function vertex_to_midpoint(fᵥ)
    return 0.5*(fᵥ[1:end-1]+fᵥ[2:end])
end