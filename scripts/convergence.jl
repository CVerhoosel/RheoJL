using RheoJL
using LaTeXStrings
using Plots, Interpolations, Statistics

# Test scenario
scenario = (τ₁ = 10.0,  # Pa
            K  = 50.0,   # Pa⋅s^n
            n  = 0.5,    # 
            η₀ = 100_000.0, # Pa⋅s
            R₀ = 0.01,   # m
            V  = 1e-6,   # m^3
            F  = 1.0,    # N
            T  = 100.0)  # s

# Observation times
tₘ =[(1/8)*2^(i-1) for i in 1:10]

# Ensure output directory exists
outdir = joinpath(@__DIR__, "..", "output")
if !isdir(outdir)
    mkpath(outdir)
end

function error(sol, sol₀, tₘ)
    itp = linear_interpolation(sol[:,1], sol[:,3])
    itp₀ = linear_interpolation(sol₀[:,1], sol₀[:,3])

    return mean( abs.( itp₀(tₘ) .- itp(tₘ) ) ./ itp₀(tₘ) )
end

function run_squeezeflow(N = 100, Δt₀ = 0.01, Nₜ = 100; scenario, tₘ)
    (; τ₁, K, n, η₀, R₀, V, F, T) = scenario
    h₀  = (V/(π*R₀^2))/2

    sol, sol_info = integrate_system(h₀, R₀, τ₁, K, n, η₀, F, N, T, Nₜ; Δt₀=Δt₀, output=:short)

    return sol, sol_info
end

function mesh_convergence(N = 2 .^(2:12), Δt₀ = 2^-9, Nₜ = 2^8, overkill=2; scenario, tₘ, T_ref)
    sols   = []
    times  = []
    
    fig_R = plot()
    for Nᵢ in N
        elapsed = @elapsed sol, sol_info = run_squeezeflow(Nᵢ, Δt₀, Nₜ; scenario=scenario, tₘ=tₘ)
        push!(times, elapsed)
        println("N=$Nᵢ: elapsed=$elapsed")  # Print elapsed time and error
        plot!(sol[:,1], sol[:,3], label="N=$Nᵢ", lw=2)
        push!(sols, sol)
    end
    display(fig_R)

    errors = []
    for sol in sols[1:end-overkill]
        push!(errors, error(sols[end], sol, tₘ))
    end

    fig_err = plot(N[1:end-overkill], errors, lw=2, xscale=:log10, yscale=:log10, xlabel=L"N_r", ylabel="Mean relative error", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    savefig(fig_err, joinpath(outdir, "error_vs_Nr.pdf"))
    display(fig_err)

    fig_time = plot(N, times/T_ref, lw=2, xscale=:log10, yscale=:log10, xlabel=L"N_t", ylabel="Normalized simulation time", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    savefig(fig_time, joinpath(outdir, "time_vs_Nr.pdf"))
    display(fig_time)
end

function time_convergence(N = 2^4, Δt₀ = 2^-9, Nₜ = 2 .^(4:14), overkill=2; scenario, tₘ, T_ref)
    sols  = []
    times = []
    
    fig_R = plot()
    for Nₜᵢ in Nₜ 
        elapsed = @elapsed sol, sol_info = run_squeezeflow(N, Δt₀, Nₜᵢ; scenario=scenario, tₘ=tₘ)
        push!(times, elapsed)
        println("Nₜ=$Nₜᵢ: elapsed=$elapsed")  # Print elapsed time and error
        plot!(sol[2:end,1], sol[2:end,3], label="Nₜ=$Nₜᵢ", lw=2)
        push!(sols, sol)
    end
    display(fig_R)

    errors = []
    for sol in sols[1:end-overkill]
        push!(errors, error(sols[end], sol, tₘ))
    end

    fig_err = plot(Nₜ[1:end-overkill], errors, lw=2, xscale=:log10, yscale=:log10, xlabel=L"N_t", ylabel="Mean relative error", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    savefig(fig_err, joinpath(outdir, "error_vs_Nt.pdf"))
    display(fig_err)

    fig_time = plot(Nₜ, times/T_ref, lw=2, xscale=:log10, yscale=:log10, xlabel=L"N_t", ylabel="Normalized simulation time", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    savefig(fig_time, joinpath(outdir, "time_vs_Nt.pdf"))
    display(fig_time)
end

function init_convergence(N = 2^4, Δt₀ = 2.0 .^(-7:-1:-16), Nₜ = 2^14, overkill=2; scenario, tₘ, T_ref)
    sols  = []
    times = []
    
    fig_R = plot()
    for Δt₀ᵢ in Δt₀
        elapsed = @elapsed sol, sol_info = run_squeezeflow(N, Δt₀ᵢ, Nₜ; scenario=scenario, tₘ=tₘ)
        push!(times, elapsed)
        println("Δt₀=$Δt₀ᵢ: elapsed=$elapsed")  # Print elapsed time and error
        plot!(sol[2:end,1], sol[2:end,3], label="Δt₀=$Δt₀ᵢ", lw=2)
        push!(sols, sol)
    end
    display(fig_R)

    errors = []
    for sol in sols[1:end-overkill]
        push!(errors, error(sols[end], sol, tₘ))
    end

    fig_err = plot(1.0 ./ Δt₀[1:end-overkill], errors, lw=2, xscale=:log10, yscale=:log10, xlabel=L"\Delta t_0^{-1}", ylabel="Mean relative error", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    savefig(fig_err, joinpath(outdir, "error_vs_Dt0.pdf"))
    display(fig_err)

    fig_time = plot(1.0 ./ Δt₀, times/T_ref, lw=2, xscale=:log10, yscale=:log10, xlabel=L"\Delta t_0^{-1}", ylabel="Normalized simulation time", legend=false, grid=:both, gridalpha=0.5, marker=:o)
    savefig(fig_time, joinpath(outdir, "time_vs_Dt0.pdf"))
    display(fig_time)
end

Nᵣ = [4, 8, 16, 32, 64, 128, 256]
Δt₀ = [2^-7, 2^-8, 2^-9, 2^-10, 2^-11, 2^-12, 2^-13]
Nₜ = [128, 256, 512, 1024, 2048, 4096, 8192]

T_ref = @elapsed ref, ref_info = run_squeezeflow(Nᵣ[end], Δt₀[end], Nₜ[end]; scenario, tₘ)
println("Reference simulation time: $T_ref [s]")

# mesh_convergence(Nᵣ, Δt₀[end], Nₜ[end], 2; scenario=scenario, tₘ=tₘ, T_ref=T_ref)
# time_convergence(Nᵣ[end], Δt₀[end], Nₜ, 2; scenario=scenario, tₘ=tₘ, T_ref=T_ref)
# init_convergence(Nᵣ[end], Δt₀, Nₜ[end], 2; scenario=scenario, tₘ=tₘ, T_ref=T_ref)

T_opt = @elapsed opt, opt_info = run_squeezeflow(16, 2^-9, 256; scenario, tₘ)
println("Optimized simulation time: $T_opt [s]")
println("Mean relative error: $(error(opt, ref, tₘ))")

TFEM = CSV.read(joinpath(@__DIR__, "..", "data", "T-FEM.txt"), DataFrame; delim=' ', ignorerepeated=true)
itp_TFEM = linear_interpolation(TFEM[!,"t[s]"], TFEM[!,"R[m]"])

fig = plot(ref[:,1], ref[:,3], label="Reference", lw=2, xlabel=L"t~[s]", ylabel=L"R~[m]")
plot!(opt[:,1], opt[:,3], label="Optimized", lw=2)
plot!(TFEM[!,"t[s]"], TFEM[!,"R[m]"], label="T-FEM", lw=2)
savefig(fig, joinpath(outdir, "R_vs_t.pdf"))
display(fig)

itp_ref = linear_interpolation(ref[:,1], ref[:,3])
itp_opt = linear_interpolation(opt[:,1], opt[:,3])
err = plot(opt[2:end,1], abs.(opt[2:end,3]-itp_ref(opt[2:end,1])) ./ itp_ref(opt[2:end,1]), label="Optimized", xscale=:log10, yscale=:log10, lw=2, color=1, xlims=(tₘ[1]/2, tₘ[end]*2), xlabel=L"t~[s]", ylabel="Relative error", legend=false, grid=:both, gridalpha=0.5, xticks=[0.1, 1, 10, 100], yticks=[0.001, 0.01], ylims=(0.001, 0.01))
plot!(tₘ, abs.(itp_opt(tₘ)-itp_ref(tₘ)) ./ itp_ref(tₘ), label="", marker=:x, markersize=5, markerstrokewidth=2, xscale=:log10, yscale=:log10, line=nothing, color=1)
savefig(err, joinpath(outdir, "error_vs_t.pdf"))
display(err)

ref_TFEM = itp_ref(TFEM[!,"t[s]"])
itp_ref_TFEM = linear_interpolation(TFEM[!,"t[s]"], ref_TFEM)
errTFEM = plot( TFEM[!,"t[s]"][2:end], abs.(TFEM[!,"R[m]"][2:end] .- ref_TFEM[2:end]) ./ TFEM[!,"R[m]"][2:end], xscale=:log10, yscale=:log10, lw=2, xlabel=L"t~[s]", ylabel="Relative error", grid=:both, gridalpha=0.5, color=1)
plot!(tₘ, abs.(itp_TFEM(tₘ)-itp_ref_TFEM(tₘ)) ./ itp_TFEM(tₘ), label="", marker=:x, markersize=5, markerstrokewidth=2, xscale=:log10, yscale=:log10, line=nothing, color=1, xlims=(tₘ[1]/2, tₘ[end]*2), legend=false)
savefig(errTFEM, joinpath(outdir, "TFEMerror_vs_t.pdf"))
display(errTFEM)