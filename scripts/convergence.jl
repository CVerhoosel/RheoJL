using RheoJL
using LaTeXStrings
using Plots, Interpolations, Statistics
using ReadVTK

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

function run_squeezeflow(N, Δt₀, Nₜ; scenario, tₘ)
    (; τ₁, K, n, η₀, R₀, V, F, T) = scenario
    h₀ = (V/(π*R₀^2))/2

    sol, sol_info = integrate_system(h₀, R₀, τ₁, K, n, η₀, F, N, T, Nₜ, Δt₀; target=nothing, Δₘₐₓ=nothing, output=:short)

    return sol, sol_info
end

function mesh_convergence(N, Δt₀, Nₜ, overkill; scenario, tₘ, T_ref)
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

function time_convergence(N, Δt₀, Nₜ, overkill=2; scenario, tₘ, T_ref)
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

function init_convergence(N, Δt₀, Nₜ, overkill; scenario, tₘ, T_ref)
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

Nᵣ = [2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8]
Δt₀ = [2^-7, 2^-8, 2^-9, 2^-10, 2^-11, 2^-12, 2^-13]
Nₜ = [2^7, 2^8, 2^9, 2^10, 2^11, 2^12, 2^13]

T_ref = @elapsed ref, ref_info = run_squeezeflow(Nᵣ[end], Δt₀[end], Nₜ[end]; scenario, tₘ)
println("Reference simulation time: $T_ref [s]")

mesh_convergence(Nᵣ, Δt₀[end], Nₜ[end], 2; scenario=scenario, tₘ=tₘ, T_ref=T_ref)
time_convergence(Nᵣ[end], Δt₀[end], Nₜ, 2; scenario=scenario, tₘ=tₘ, T_ref=T_ref)
init_convergence(Nᵣ[end], Δt₀, Nₜ[end], 2; scenario=scenario, tₘ=tₘ, T_ref=T_ref)

T_opt = @elapsed opt, opt_info = run_squeezeflow(Nᵣ[2], Δt₀[2], Nₜ[3]; scenario, tₘ)
println("Optimized simulation time: $T_opt [s]")
println("Mean relative error: $(error(opt, ref, tₘ))")

TFEM = load_data("T-FEM.txt")
itp_TFEM = linear_interpolation(TFEM[!,"t[s]"], TFEM[!,"R[m]"])

fig = plot(ref[:,1], ref[:,3], label="Reference", lw=2, xlabel=L"t~[s]", ylabel=L"R~[m]")
plot!(opt[:,1], opt[:,3], label="Optimized", lw=2)
plot!(TFEM[!,"t[s]"], TFEM[!,"R[m]"], label="T-FEM", lw=2)
savefig(fig, joinpath(outdir, "R_vs_t.pdf"))
display(fig)

itp_ref = linear_interpolation(ref[:,1], ref[:,3])
itp_opt = linear_interpolation(opt[:,1], opt[:,3])
err = plot(opt[2:end,1], abs.(opt[2:end,3]-itp_ref(opt[2:end,1])) ./ itp_ref(opt[2:end,1]), label="Optimized", xscale=:log10, yscale=:log10, lw=2, color=1, xlims=(tₘ[1]/2, tₘ[end]*2), xlabel=L"t~[s]", ylabel="Relative error", legend=false, grid=:both, gridalpha=0.5, xticks=[0.1, 1, 10, 100], yticks=[0.001, 0.01], ylims=(0.001, 0.1))
plot!(tₘ, abs.(itp_opt(tₘ)-itp_ref(tₘ)) ./ itp_ref(tₘ), label="", marker=:x, markersize=5, markerstrokewidth=2, xscale=:log10, yscale=:log10, line=nothing, color=1)
savefig(err, joinpath(outdir, "error_vs_t.pdf"))
display(err)

ref_TFEM = itp_ref(TFEM[!,"t[s]"])
itp_ref_TFEM = linear_interpolation(TFEM[!,"t[s]"], ref_TFEM)
errTFEM = plot( TFEM[!,"t[s]"][2:end], abs.(TFEM[!,"R[m]"][2:end] .- ref_TFEM[2:end]) ./ TFEM[!,"R[m]"][2:end], xscale=:log10, yscale=:log10, lw=2, xlabel=L"t~[s]", ylabel="Relative error", grid=:both, gridalpha=0.5, color=1)
plot!(tₘ, abs.(itp_TFEM(tₘ)-itp_ref_TFEM(tₘ)) ./ itp_TFEM(tₘ), label="", marker=:x, markersize=5, markerstrokewidth=2, xscale=:log10, yscale=:log10, line=nothing, color=1, xlims=(tₘ[1]/2, tₘ[end]*2), legend=false)
savefig(errTFEM, joinpath(outdir, "TFEMerror_vs_t.pdf"))
display(errTFEM)

tᵢ = opt[end,1]
Rᵢ = opt[end,3]
hᵢ = scenario.V/(2*π*Rᵢ^2)
∂h∂tᵢ, infoᵢ = solve_system(hᵢ, Rᵢ, scenario.τ₁, scenario.K, scenario.n, scenario.η₀, scenario.F, Nᵣ[3]; output=:long)
rᵢ = vertex_to_midpoint(collect(range(0, Rᵢ, length=Nᵣ[3])))
Qᵢ = -∂h∂tᵢ * rᵢ
∂p∂rᵢ = first.(solve_∂p∂r.(hᵢ, scenario.τ₁, scenario.K, scenario.n, scenario.η₀, Qᵢ; bracket=[-Inf,0]))
pᵢ = right_integrate(collect(range(0, Rᵢ, length=Nᵣ[3])), ∂p∂rᵢ)

∂p∂r_interp = linear_interpolation(rᵢ, ∂p∂rᵢ)

#########################
# Comparison with T-FEM #
#########################

# Load VTK file
h5_file = joinpath(@__DIR__, "..", "data", "T-FEM-t100.h5")

using JLD2
data = jldopen(h5_file)
points = data["points"][1:2, :]'

nr = 41
nz = 41

r = points[1:nr:end,2]
z = points[1:nz,1]
z = z .- z[1]
v = zeros(nz, nr)
τ = zeros(nz, nr)
for ir in 1:nr
   for iz in 1:nz
         idx = (ir-1)*nz + nz + 1 - iz
         v[iz,ir] = data["point_data"]["velocity"][2,idx]
         τ[iz,ir] = data["point_data"]["tau_zr"][idx]
    end
end

vc_fig = contourf(r ./ Rᵢ, z ./ hᵢ, v, 
    levels=16,
    color=:turbo, 
    xlabel=L"r~[mm]", 
    ylabel=L"z~[mm]",
    colorbar_title=L"v~[mm/s]",
    colorbar=:bottom,
    colorbar_formatter=:scientific,
    fillalpha=1,
    linewidth=0.0
    )

savefig(vc_fig, joinpath(outdir, "vcontour.pdf"))
display(vc_fig)

τc_fig = contourf(r ./ Rᵢ, z ./ hᵢ, τ / scenario.τ₁, 
    levels=collect(0:0.25:3.5),
    color=:turbo, 
    xlabel=L"r / R", 
    ylabel=L"2 z / H",
    colorbar_title=L"τ / \tau_y",
    colorbar=:bottom,
    colorbar_formatter=:scientific,
    fillalpha=1.0,
    linewidth=0.0
    )

plot!(rᵢ / Rᵢ, (scenario.τ₁ ./ abs.(∂p∂rᵢ)) / hᵢ, lw=4, color=:white, label="", xlims=[0,1], ylims=[0,1],
    marker=:circle, markersize=5, markerstrokewidth=0)

savefig(τc_fig, joinpath(outdir, "τcontour.pdf"))
display(τc_fig)

rslices = [11, 21, 31]


v_plot = plot()
for (i, rslice) in enumerate(rslices)
    plot!(v[:,rslice] ./ v[1,end], z ./ hᵢ, lw=2, ylabel=L"2 z / H", xlabel=L"v / v(R,0)", label="r / R = $(round(r[rslice]/Rᵢ, digits=2))", color=i, ls=:dash)

    zs, vs, np1 = velocity_profile(hᵢ, ∂p∂r_interp(r[rslice]), scenario.τ₁, scenario.K, scenario.n, scenario.η₀)
    plot!(vs ./ v[1,end], zs ./ hᵢ, lw=2, label="", color=i)
end
savefig(v_plot, joinpath(outdir, "vprofiles.pdf"))
display(v_plot)