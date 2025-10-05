using RheoJL
using LaTeXStrings
using Plots

# Test scenario
scenario = (τ₁ = 100.0,  # Pa
            K  = 50.0,   # Pa⋅s^n
            n  = 0.5,    # 
            η₀ = 1000.0, # Pa⋅s
            R₀ = 0.01,   # m
            V  = 1e-6,   # m^3
            F  = 1.0,    # N
            T  = 100.0)  # s

# Observation times
t =[(1/8)*2^(i-1) for i in 1:10]

function run_squeezeflow(N = 100, Δ₀ = 0.01, Δtₘₐₓ = 1; scenario)
    (; τ₁, K, n, η₀, R₀, V, F, T) = scenario
    h₀  = (V/(π*R₀^2))/2

    sol, sol_info = integrate_system(h₀, R₀, τ₁, K, n, η₀, F, N, T; Δ₀=Δ₀, Δtₘₐₓ=Δtₘₐₓ, output=:long)

    plot(sol[:,1], sol[:,3], label="", lw=3, xlabel=L"t~[s]", ylabel=L"R~[m]")

    # display(plt)
    # # Single time step
    # ∂h∂t, info = solve_system(h₀, R₀, τ₁, K, n, η₀, F, N; output=:long)

    # for infoᵢ in eachrow(info)
    #     println(infoᵢ)
    # end
end

run_squeezeflow(; scenario=scenario)