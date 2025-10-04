using RheoJL

function main(
    τ₁ = 100.0,
    K  = 1.0,
    n  = 1.0,
    η₀ = 1000.0,
    R₀ = 0.013,
    V  = 0.73e-6,
    F  = 0.5,
    N  = 100)

    h₀  = (V/(π*R₀^2))/2

    # Single time step
    ∂h∂t, info = solve_system(h₀, R₀, τ₁, K, n, η₀, F, N; output=:long, outputinner=:short)

    for infoᵢ in eachrow(info)
        println(infoᵢ)
    end
end