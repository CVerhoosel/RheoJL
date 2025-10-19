using RheoJL
using Distributions, Random, Statistics
using ProgressMeter, StatsPlots, LaTeXStrings

n_lower = 0.1
n_upper = 1.0
nbins = 9

# Test scenario
scenario = (τ₁ = Uniform(0.0, 100.0),  # Pa
            K  = truncated(Normal(50.0, 25.0); lower=0.0),   # Pa⋅s^n
            n  = truncated(Normal(0.5, 0.25);lower=n_lower, upper=n_upper), 
            log₁₀η₀ = Uniform(2, 7), # Pa⋅s
            R₀ = truncated(Normal(1e-2, 0.5e-2);lower=0.001), # m
            V  = truncated(Normal(1e-6, 0.5e-6);lower=0.001), # m^3
            F  = truncated(Normal(1.0, 0.5);lower=0.0), # N
            T  = 100.0) # s

function scenario_mean(scenario::NamedTuple)
    names = Tuple(keys(scenario))
    vals = tuple( ((v isa Distribution ? mean(v) : v) for v in values(scenario))... )
    return NamedTuple{names}(vals)
end

function scenario_rand(scenario::NamedTuple)
    names = Tuple(keys(scenario))
    vals = tuple( ((v isa Distribution ? rand(v) : v) for v in values(scenario))... )
    return NamedTuple{names}(vals)
end

function run_squeezeflow(Nᵣ, Δt₀, Nₜ; scenario)
    (; τ₁, K, n, log₁₀η₀, R₀, V, F, T) = scenario
    h₀ = (V/(π*R₀^2))/2
    η₀ = 10.0^log₁₀η₀

    sol, sol_info = integrate_system(h₀, R₀, τ₁, K, n, η₀, F, Nᵣ, T, Nₜ, Δt₀; output=:long)

    info = [reduce(vcat,[iter_info[5:end] for iter_info in sol_infoᵢ]) for sol_infoᵢ in sol_info]
    
    return sol, [τ₁, K, n, log₁₀η₀, R₀, V, F, sum(size.(info,1))]
end

Nₛ = 10_000  # Number of samples
Nᵣ = 2^4 # Number of radial points
Nₜ = 2^9 # Number of time steps
Δt₀ = 2.0^-9 # Initial time step

sol₀, info₀ = run_squeezeflow(Nᵣ, Δt₀, Nₜ; scenario=scenario_mean(scenario))
neval₀ = info₀[end]
sample = []

@showprogress for i in 1:Nₛ
    scenarioᵢ = scenario_rand(scenario)
    sol, info = run_squeezeflow(Nᵣ, Δt₀, Nₜ; scenario=scenarioᵢ)
    info[end] /= neval₀
    push!(sample, info)
end

sample_mat = reduce(vcat,sample')

function bin_iterations_by_col(sample_mat::AbstractMatrix, col::Integer, edges)

    nbins = length(edges) - 1

    # prepare bins: each bin will collect the iterations (last column) as Float64
    bins = [Float64[] for _ in 1:nbins]
    for r in 1:size(sample_mat,1)
        v = sample_mat[r, col]
        idx = searchsortedlast(edges, v)
        push!(bins[idx], sample_mat[r, end])
    end

    return bins
end

edges = collect(range(n_lower, n_upper; length=nbins+1))
centers = 0.5*(edges[2:end]+edges[1:end-1])
binned = bin_iterations_by_col(sample_mat, 3, edges)

# filter empty bins first
nonempty = findall(!isempty, binned)
binned_nonempty = binned[nonempty]
centers_nonempty = centers[nonempty]

# category indices 1..k
idxs = 1:length(binned_nonempty)

boxplot(idxs, binned_nonempty;
    width = 0.5,
    whisker_width = :match,
    xlabel = "n",
    ylabel = "Normalized number of function evaluations",
    legend = false,
    grid = :both,
    gridalpha = 0.5,
    xticks = (idxs, string.(round.(centers_nonempty, digits=3))))