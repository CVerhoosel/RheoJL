using RheoJL
using LaTeXStrings
using Plots, Interpolations, Statistics

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
tₘ =[(1/8)*2^(i-1) for i in 1:10]

function error(sol, sol₀, tₘ)
    itp = linear_interpolation(sol[:,1], sol[:,3])
    itp₀ = linear_interpolation(sol₀[:,1], sol₀[:,3])

    return mean( abs.( itp₀(tₘ) .- itp(tₘ) ) ./ itp₀(tₘ) )
end

function run_squeezeflow(N = 100, Δt₀ = 0.01, Δtₘₐₓ = 1; scenario, tₘ)
    (; τ₁, K, n, η₀, R₀, V, F, T) = scenario
    h₀  = (V/(π*R₀^2))/2

    sol, sol_info = integrate_system(h₀, R₀, τ₁, K, n, η₀, F, N, T; Δt₀=Δt₀, Δtₘₐₓ=Δtₘₐₓ, output=:long)

    return sol, sol_info
end

function mesh_convergence(N = [4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048], Δt₀ = 1e-6, Δtₘₐₓ = 1, overkill=3; scenario, tₘ)
    sols   = []
    times  = []
    
    fig_R = plot()
    for Nᵢ in N
        elapsed = @elapsed sol, sol_info = run_squeezeflow(Nᵢ, Δt₀, Δtₘₐₓ; scenario=scenario, tₘ=tₘ)
        push!(times, elapsed)
        println("N=$Nᵢ: elapsed=$elapsed")  # Print elapsed time and error
        plot!(sol[:,1], sol[:,3], label="N=$Nᵢ", lw=2)
        push!(sols, sol)
    end
    display(fig_R)

    errors = []
    for (i, sol) in enumerate(sols[1:end-overkill])
        push!(errors, error(sols[end], sol, tₘ))
    end

    fig_err = plot(N[1:end-overkill], errors, lw=2, xscale=:log10, yscale=:log10, xlabel=L"N", ylabel=L"e_R", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    display(fig_err)

    fig_time = plot(N, times / times[end], lw=2, xscale=:log10, yscale=:log10, xlabel=L"N", ylabel=L"T_{\rm sim}", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    display(fig_time)
end

function init_convergence(N = 128, Δt₀ = [2^-7, 2^-8, 2^-9, 2^-10, 2^-11, 2^-12, 2^-13, 2^-14, 2^-15], Δtₘₐₓ = 100., overkill=4; scenario, tₘ)
    sols  = []
    times = []
    steps = []
    
    fig_R = plot()
    for Δt₀ᵢ in Δt₀
        elapsed = @elapsed sol, sol_info = run_squeezeflow(N, Δt₀ᵢ, Δtₘₐₓ; scenario=scenario, tₘ=tₘ)
        push!(times, elapsed)
        println("Δt₀=$Δt₀ᵢ: elapsed=$elapsed")  # Print elapsed time and error
        plot!(sol[:,1], sol[:,3], label="Δt₀=$Δt₀ᵢ", lw=2)
        push!(sols, sol)
        push!(steps, size(sol,1))
    end
    display(fig_R)

    errors = []
    fig_errR = plot(xscale=:log10)
    for (i, sol) in enumerate(sols[1:end-overkill])

        itp = linear_interpolation(sol[:,1], sol[:,3])
        itp₀ = linear_interpolation(sols[end][:,1], sols[end][:,3])
        plot!(tₘ, (itp(tₘ)-itp₀(tₘ))./itp₀(tₘ), label="Δt₀=$(Δt₀[i])", color=i)


        for j in 1:size(sol,1)-1
            if abs( (sol[j+1,1]-sol[j,1]) - Δtₘₐₓ ) < 1e-6*Δtₘₐₓ
                plot!([sol[j,1]], [(sol[j,3]-itp₀(sol[j,1]))/itp₀(sol[j,1])], color=i, marker=:o, label="")
                break
            end
        end

        push!(errors, error(sols[end], sol, tₘ))
    end
    display(fig_errR)

    fig_err = plot(1 ./ Δt₀[1:end-overkill], errors, lw=2, xscale=:log10, yscale=:log10, xlabel=L"\Delta t_0^{-1}", ylabel=L"e_R", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    display(fig_err)

    fig_time = plot(1 ./ Δt₀, times / times[end], lw=2, xscale=:log10, yscale=:log10, xlabel=L"\Delta t_0^{-1}", ylabel=L"T_{\rm sim}", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    display(fig_time)
end

function time_convergence(N = 128, Δt₀ = 1e-6, Δtₘₐₓ = [2^1, 2^0, 2^-1, 2^-2,2^-3, 2^-4, 2^-5, 2^-6, 2^-7, 2^-8], overkill=1; scenario, tₘ)
    sols  = []
    times = []
    steps = []
    
    fig_R = plot(xscale=:log10, yscale=:log10)
    for Δtₘₐₓᵢ in Δtₘₐₓ
        elapsed = @elapsed sol, sol_info = run_squeezeflow(N, Δt₀, Δtₘₐₓᵢ; scenario=scenario, tₘ=tₘ)
        push!(times, elapsed)
        println("Δtₘₐₓ=$Δtₘₐₓᵢ: elapsed=$elapsed")  # Print elapsed time and error
        plot!(sol[2:end,1], sol[2:end,3], label="Δtₘₐₓ=$Δtₘₐₓᵢ", lw=2)
        push!(sols, sol)
        push!(steps, size(sol,1))
    end
    display(fig_R)

    errors = []
    fig_errR = plot(xscale=:log10)
    for (i, sol) in enumerate(sols[1:end-overkill])

        itp = linear_interpolation(sol[:,1], sol[:,3])
        itp₀ = linear_interpolation(sols[end][:,1], sols[end][:,3])
        plot!(tₘ, (itp(tₘ)-itp₀(tₘ))./itp₀(tₘ), label="Δtₘₐₓ=$(Δtₘₐₓ[i])", color=i)


        for j in 1:size(sol,1)-1
            if abs( (sol[j+1,1]-sol[j,1]) - Δtₘₐₓ[i] ) < 1e-6*Δtₘₐₓ[i]
                plot!([sol[j,1]], [(sol[j,3]-itp₀(sol[j,1]))/itp₀(sol[j,1])], color=i, marker=:o, label="")
                break
            end
        end

        push!(errors, error(sols[end], sol, tₘ))
    end
    display(fig_errR)

    fig_err = plot(1 ./ Δtₘₐₓ[1:end-overkill], errors, lw=2, xscale=:log10, yscale=:log10, xlabel=L"\Delta t_{\rm max}^{-1}", ylabel=L"e_R", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    display(fig_err)

    fig_time = plot(1 ./ Δtₘₐₓ, times / times[end], lw=2, xscale=:log10, yscale=:log10, xlabel=L"\Delta t_{\rm max}^{-1}", ylabel=L"T_{\rm sim}", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    display(fig_time)
end

# mesh_convergence(; scenario=scenario, tₘ=tₘ)
time_convergence(; scenario=scenario, tₘ=tₘ)