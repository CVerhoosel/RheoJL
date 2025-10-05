### A Pluto.jl notebook ###
# v0.20.16

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 791fbd6c-1a96-4919-ac04-6ee5abc1a538
# ╠═╡ show_logs = false
begin
	using Pkg
	
    # activate the global environment
    Pkg.activate()
	
	using Revise, PlutoUI, Plots, LaTeXStrings, Statistics, StatsPlots
	using RheoJL

md"""
!!! info "Julia environment"
	This notebook makes use of the global environment. Make sure all used packages, including RheoJL, are available there.
"""
end

# ╔═╡ 316a838d-57db-4194-b4e1-857f3afc44f6
md"""
# RheoJL single time step test
"""

# ╔═╡ 2d4046cd-fd85-4f02-8126-469d296fda1a
md"## Parameters"

# ╔═╡ 8c337acd-42f9-4850-8f6f-d3397a3183c1
@bind reset_parameters CounterButton("Reset parameters")

# ╔═╡ 75987388-9b29-4505-b641-e1d4c8d72918
md"### Fluid"

# ╔═╡ 8c15239f-c3e8-4e79-99a8-9b8c5b56bc25
begin
	reset_parameters
	md"""
	``\tau_{y}~[Pa]``: $(@bind τ₁ Slider(0:1:300; default=130, show_value=true))
	"""
end

# ╔═╡ c58882b9-ca14-479f-ad4d-a98e194b5774
begin
	reset_parameters
	md"""
	``K~[Pa\cdot s^n]``: $(@bind K Slider(1:100; default=30, show_value=true))
	"""
end

# ╔═╡ a8ba3e65-e0a2-4a9f-bda4-5c4a4a64f7a1
begin
	reset_parameters
	md"""
	``n~[-]``: $(@bind n Slider(0.01:0.01:1.5; default=0.4, show_value=true))
	"""
end

# ╔═╡ 8c36449a-99f8-4bbc-9194-c432aae6708b
begin
	reset_parameters
	md"""
	``\eta_0~[Pa \cdot s]``: $(@bind η₀ Slider(10:10:10_000; default=1000, show_value=true))
	"""
end

# ╔═╡ 99a62371-6939-4c1f-af94-eeb7d177ceb7
begin
	#ε = 10.0^ε_input
	#η₀ = 10*((n*K*τ₁^(n-1))/ε)^(1/n)
	∂γ∂t₀ = τ₁/η₀
	τ₀ = τ₁ - K*∂γ∂t₀^n
	md"""
	| Parameter | Value | Unit |
	|-----------|-------|------|
	| ``\tau_y`` | $(τ₁) | ``Pa`` |
	| ``\tau_0`` | $(τ₀) | ``Pa`` |
	| ``K`` | $(K) | ``Pa \cdot s^n`` |
	| ``n`` | $(n) | ``-`` |
	| ``\eta_0`` | $(η₀) | ``Pa \cdot s`` |
	"""
end

# ╔═╡ f3407c7b-4686-4a9d-b354-5d3e9c8a4033
md"### System"

# ╔═╡ dab6942e-3ffa-458c-b09c-a90dd2cd622b
md"""
Allow slip: $(@bind allow_slip CheckBox(default=false))
"""

# ╔═╡ ef45baff-f705-488b-a0e5-f53c5351f452
md"""
Include capillary pressure: $(@bind include_capillary CheckBox(default=false))
"""

# ╔═╡ 3a3ead59-460a-4730-a012-c1244d163446
begin
	reset_parameters
	md"``R_0~[mm]``: $(@bind R₀_input Slider(5:0.1:25; default=13, show_value=true))"
end

# ╔═╡ 584713b9-455a-4df4-9694-de8acf84f801
begin
	reset_parameters
	md"""
	``V~[ml]``: $(@bind V_input Slider(0.1:0.01:1.5; default=0.70, show_value=true))
	"""
end

# ╔═╡ 25f30cff-5b1a-4309-99b2-1e59a68d69c5
begin
	reset_parameters
	md"""
	``F~[N]``: $(@bind F Slider(0.1:0.1:5.0; default=1.5, show_value=true))
	"""
end

# ╔═╡ b87384d8-97c1-4c55-a43d-0e7c16394f2b
begin
	if allow_slip
	reset_parameters
	md"""
	``\log_{10} \beta~[Pa \cdot s /m]``: $(@bind β_input Slider(0:0.2:10; default=6, show_value=true))
	"""
	end
end

# ╔═╡ e1c4d2bd-342b-41ad-a2ef-75e50614140d
begin
	if include_capillary
	reset_parameters
	md"""
	``\gamma~[N/m]``: $(@bind γ_input Slider(0:0.01:0.3; default=0.06, show_value=true))
	"""
	end
end

# ╔═╡ a6f0d7c4-40ea-4d02-9bdc-ff5d1d016d49
begin
	if include_capillary
	reset_parameters
	md"""
	``\alpha~[-]``: $(@bind α_input Slider(0:0.1:1; default=0.5, show_value=true))
	"""
	end
end

# ╔═╡ 2d129381-d823-403f-8eec-4322da51ddd7
begin
	V  = V_input / 1e6
	R₀ = R₀_input / 1e3
	h  = (V/(π*R₀^2))/2
	β  = allow_slip ? 10.0^β_input : nothing
	γ  = include_capillary ? γ_input : nothing
	α  = include_capillary ? α_input : nothing
	md"""
	| Parameter | Value | Unit |
	|-----------|-------|------|
	| ``R_0`` | $(R₀) | ``m`` |
	| ``V``   | $(V) | ``m^3`` |
	| ``h_0`` | $(h) | ``m`` |
	| ``F``   | $(F) | ``N`` |
	| ``β``   | $(β) | ``Pa\cdot s / m`` |
	| ``γ``   | $(γ) | ``N / m`` |
	| ``α``   | $(α) | ``-`` |
	"""
end

# ╔═╡ c0a98c3e-4f48-4d3a-86de-e769a7fdf502
md"### Discretization"

# ╔═╡ da4031fa-1ff6-45ff-ab48-2ba31a6f1300
begin
	reset_parameters
	md"""
	``N``: $(@bind N Slider(10:1000; default=100, show_value=true))
	"""
end

# ╔═╡ 76aa1611-f856-4baa-b6c5-095d85c86445
begin
	reset_parameters
	md"""
	``T~[s]``: $(@bind T Slider(10:10:1_000; default=150, show_value=true))
	"""
end

# ╔═╡ 9cc0c580-652e-44d8-8e0d-7291fc16e36e
begin
	reset_parameters
	md"""
	``\Delta_0``: $(@bind Δ₀ Slider(0.001:0.001:0.1; default=0.01, show_value=true))
	"""
end

# ╔═╡ ffe187b0-45e6-4b42-98c8-5013d232ee01
begin
	rᵥ = collect(range(0, R₀, length=N))
    rₘ = vertex_to_midpoint(rᵥ)
	md"""
	| Parameter | Value | Unit |
	|-----------|-------|------|
	| ``N`` | $(N) | ``-`` |
	| ``T`` | $(T) | ``s`` |
	| ``\Delta_0=-\frac{\dot{h}_0 \Delta t_0 }{ h_0 }`` |   $(Δ₀) | ``-`` |
	"""
end

# ╔═╡ 41529b16-8175-4c8d-9f0f-8ba63f83ff6b
begin
	reset_parameters
	md"""
	``\Delta t_{\rm max}``: $(@bind Δtₘₐₓ Slider(1:1:100; default=10, show_value=true))
	"""
end

# ╔═╡ e024bd4c-e501-42c8-9396-6900b3f5c583
md"""## Time step solver"""

# ╔═╡ f4fe528d-d827-49fb-a6a9-59c5508aadcb
begin
	∂h∂t, info = solve_system(h, R₀, τ₁, K, n, η₀, F, N; β=β, γ=γ, α=α, output=:long)
	∂h∂t_range = collect(range(2*∂h∂t, 0, length=100))
	F_range = []
	for ∂h∂t ∈ ∂h∂t_range
		Qₘ = -∂h∂t*rₘ
		∂p∂r₀ₘ = -(4*F)/(π*R₀^4) * rₘ
		results = [solve_∂p∂r(h, τ₁, K, n, η₀, Qₘ[i]; β=β, ∂p∂r₀=∂p∂r₀ₘ[i]) for i in 	eachindex(Qₘ)]
		∂p∂rₘ = first.(results)
		push!(F_range, force(rᵥ, ∂p∂rₘ; h=h, γ=γ, α=α))
	end
	plot(∂h∂t_range, F_range, label="", xlabel=L"\dot{h}", ylabel=L"F")
	plot!([∂h∂t_range[1], ∂h∂t_range[end]], [F, F], label="Target")
	plot!([∂h∂t],[F], label="Solution", m=:star, markersize=7)

	for (i, (∂h∂tᵢ, Fᵢ, left, right, iterations...)) in enumerate(info)
		plot!([∂h∂tᵢ], [Fᵢ], m=:circle, label="Iteration $(i)", markersize=3)
	end
	plot!()
end

# ╔═╡ 3ab29162-283d-44d0-a378-0f93937e3f3f
	boxplot(reduce(hcat,info)[5:end,:], label="", xlabel=L"i_{\rm outer}", ylabel=L"n_{\rm inner}")

# ╔═╡ f339ac8d-4544-41fd-ae00-7cac56c49215
md"""## Time integration"""

# ╔═╡ 7ffd160c-3125-41ff-bb6c-b72fa26f4e08
@bind data_files MultiCheckBox(list_data_files())

# ╔═╡ 8a353e56-f77f-4f9d-b2eb-26a04507db7c
begin
	sol, sol_info = integrate_system(h, R₀, τ₁, K, n, η₀, F, N, T; β=β, γ=γ, α=α, Δ₀=Δ₀, Δtₘₐₓ=Δtₘₐₓ, output=:long, progress=true, targetiter=5)
	plot(sol[:,1], sol[:,3], label="Model", lw=3, xlabel=L"t~[s]", ylabel=L"R~[m]")

	for data_file in data_files
		df = load_data(data_file)

		t = Vector(df[:,:Time_1])
		
		radius_cols = filter(name -> startswith(String(name), "Radius"), names(df))
		df_radius = df[:, radius_cols]

		μ = mean(Matrix(df_radius), dims=2)
		if length(radius_cols) > 1
			σ = std(Matrix(df_radius), dims=2)
			plot!(t, μ, yerr=σ, label=data_file)
		else
			plot!(t, μ, label=data_file)
		end
	end
	plot!()
end

# ╔═╡ 9faeebaf-4210-402e-b903-2b1f2e2ed5a9
let
	bar( [size(sol_infoᵢ,1)-1 for sol_infoᵢ in sol_info], label="", xlabel=L"i_{\rm time}", ylabel=L"n_{\rm outer}" )
	plot!(twinx(), sol[2:end,1]-sol[1:end-1,1], color=:red, ylabel=L"\Delta t~[s]", yscale=:log10, label="", linewidth=2)
end

# ╔═╡ 27577ae2-a446-49a8-a8f1-68198d89718d
let
	red_info = [reduce(vcat,[iter_info[5:end] for iter_info in sol_infoᵢ]) for sol_infoᵢ in sol_info]
	boxplot(red_info, label="", xlabel=L"i_{\rm time}", ylabel=L"n_{\rm inner}")
end

# ╔═╡ Cell order:
# ╟─316a838d-57db-4194-b4e1-857f3afc44f6
# ╟─791fbd6c-1a96-4919-ac04-6ee5abc1a538
# ╟─2d4046cd-fd85-4f02-8126-469d296fda1a
# ╟─8c337acd-42f9-4850-8f6f-d3397a3183c1
# ╟─75987388-9b29-4505-b641-e1d4c8d72918
# ╟─99a62371-6939-4c1f-af94-eeb7d177ceb7
# ╟─8c15239f-c3e8-4e79-99a8-9b8c5b56bc25
# ╟─c58882b9-ca14-479f-ad4d-a98e194b5774
# ╟─a8ba3e65-e0a2-4a9f-bda4-5c4a4a64f7a1
# ╟─8c36449a-99f8-4bbc-9194-c432aae6708b
# ╟─f3407c7b-4686-4a9d-b354-5d3e9c8a4033
# ╟─2d129381-d823-403f-8eec-4322da51ddd7
# ╟─dab6942e-3ffa-458c-b09c-a90dd2cd622b
# ╟─ef45baff-f705-488b-a0e5-f53c5351f452
# ╟─3a3ead59-460a-4730-a012-c1244d163446
# ╟─584713b9-455a-4df4-9694-de8acf84f801
# ╟─25f30cff-5b1a-4309-99b2-1e59a68d69c5
# ╟─b87384d8-97c1-4c55-a43d-0e7c16394f2b
# ╟─e1c4d2bd-342b-41ad-a2ef-75e50614140d
# ╟─a6f0d7c4-40ea-4d02-9bdc-ff5d1d016d49
# ╟─c0a98c3e-4f48-4d3a-86de-e769a7fdf502
# ╟─ffe187b0-45e6-4b42-98c8-5013d232ee01
# ╟─da4031fa-1ff6-45ff-ab48-2ba31a6f1300
# ╟─76aa1611-f856-4baa-b6c5-095d85c86445
# ╟─9cc0c580-652e-44d8-8e0d-7291fc16e36e
# ╟─41529b16-8175-4c8d-9f0f-8ba63f83ff6b
# ╟─e024bd4c-e501-42c8-9396-6900b3f5c583
# ╟─f4fe528d-d827-49fb-a6a9-59c5508aadcb
# ╟─3ab29162-283d-44d0-a378-0f93937e3f3f
# ╟─f339ac8d-4544-41fd-ae00-7cac56c49215
# ╟─7ffd160c-3125-41ff-bb6c-b72fa26f4e08
# ╠═8a353e56-f77f-4f9d-b2eb-26a04507db7c
# ╟─9faeebaf-4210-402e-b903-2b1f2e2ed5a9
# ╟─27577ae2-a446-49a8-a8f1-68198d89718d
