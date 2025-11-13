### A Pluto.jl notebook ###
# v0.20.21

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
	
	using Revise, PlutoUI, Plots, LaTeXStrings
	using RheoJL

md"""
!!! info "Julia environment"
	This notebook makes use of the global environment. Make sure all used packages, including RheoJL, are available there.
"""
end

# ╔═╡ 316a838d-57db-4194-b4e1-857f3afc44f6
md"""
# RheoJL rheology notebook
"""

# ╔═╡ 2d4046cd-fd85-4f02-8126-469d296fda1a
md"## Parameters"

# ╔═╡ 8c337acd-42f9-4850-8f6f-d3397a3183c1
@bind reset_parameters CounterButton("Reset parameters")

# ╔═╡ 75987388-9b29-4505-b641-e1d4c8d72918
md"### Fluid"

# ╔═╡ dd4aebde-8087-402b-b17f-59b17522cda8
begin
	reset_parameters
	models = ["Newton (K-variant)", "Newton (η₀-variant)", "Bingham", "Herschel-Bulkley", "Biviscous", "Biviscous power law"]
	@bind model Select(models, default=models[end])
end

# ╔═╡ f3407c7b-4686-4a9d-b354-5d3e9c8a4033
md"### System"

# ╔═╡ 8c15239f-c3e8-4e79-99a8-9b8c5b56bc25
begin
	reset_parameters
	md"""
	``\tau_{y}~[Pa]``: $(@bind input_τ₁ Slider(0:1:100; default=10, show_value=true))
	"""
end

# ╔═╡ c58882b9-ca14-479f-ad4d-a98e194b5774
begin
	reset_parameters
	md"""
	``K~[Pa\cdot s^n]``: $(@bind input_K Slider(1:100; default=50, show_value=true))
	"""
end

# ╔═╡ a8ba3e65-e0a2-4a9f-bda4-5c4a4a64f7a1
begin
	reset_parameters
	md"""
	``n~[-]``: $(@bind input_n Slider(0.1:0.01:1.5; default=0.5, show_value=true))
	"""
end

# ╔═╡ 8c36449a-99f8-4bbc-9194-c432aae6708b
begin
	reset_parameters
	md"""
	``{\rm log}_{10} \eta_0~[Pa \cdot s]``: $(@bind logη₀ Slider(1:8; default=5, show_value=true))
	"""
end

# ╔═╡ 99a62371-6939-4c1f-af94-eeb7d177ceb7
begin

	# K
	if model ∈ models[[1, 3, 4, 5, 6]]
		K  = input_K
	else
		K  = nothing
	end
	
	# Yield stress & K
	if model ∈ models[[3, 4, 5, 6]]
		τ₁ = input_τ₁
	else
		τ₁ = nothing
	end

	# Power law index
	if model ∈ models[[4, 6]]
		n = input_n
	else
		n = nothing
	end
	
	# Pre-yield viscosity
	if model ∈ models[[2, 5, 6]]
		η₀ = 10^logη₀
	else
		η₀ = nothing	
	end

	if K!==nothing && τ₁!==nothing && η₀!==nothing && n!==nothing
		τ₀ = τ₁ - K*(τ₁/η₀)^n
	else
		τ₀ = nothing
	end
	
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

# ╔═╡ 256361ad-327e-448c-ba72-cafccd4cea77
begin
	reset_parameters
	md"""
	``\log_{10} \dot{\gamma}_{\rm min}~[1/s]``: $(@bind log_γmin Slider(-5:0.2:1; default=-4, show_value=true))
	"""
end

# ╔═╡ 2bf18333-273f-493c-93ae-0e0a9ae4eefe
begin
	reset_parameters
	md"""
	``\log_{10} \dot{\gamma}_{\rm max}~[1/s]``: $(@bind log_γmax Slider(log_γmin+1:0.2:log_γmin+4; default=log_γmin+3, show_value=true))
	"""
end

# ╔═╡ e37fa6c6-0fb5-41ec-9818-632f245131f7
begin
	reset_parameters
	md"""
	``N_{\dot{\gamma}}``: $(@bind Nγ Slider(10:10:1_000; default=100, show_value=true))
	"""
end

# ╔═╡ 3772106b-1e2f-41a6-94bd-64a9459239c2
md"""
Allow slip: $(@bind allow_slip CheckBox(default=true))
"""

# ╔═╡ 3a3ead59-460a-4730-a012-c1244d163446
begin
	reset_parameters
	if allow_slip
		md"``H~[mm]``: $(@bind input_H Slider(5:0.1:25; default=10, show_value=true))"
	end
end

# ╔═╡ 0fc3ad2f-6742-4b08-9114-b8d831498463
begin
	reset_parameters
	if allow_slip
	md"""
	``\log_{10} \beta~[Pa \cdot s /m]``: $(@bind input_β Slider(0:0.2:10; default=6, show_value=true))
	"""
	end
end

# ╔═╡ 2d129381-d823-403f-8eec-4322da51ddd7
begin
	H = allow_slip ? 1e-3*input_H : nothing
	β = allow_slip ? 10.0^input_β : nothing
	γmin = 10.0^log_γmin
	γmax = 10.0^log_γmax
	γrng = geometric_time_sequence(γmax, γmin, Nγ)[2:end]
	md"""
	| Parameter | Value | Unit |
	|-----------|-------|------|
	| ``\dot{\gamma}_{\rm min}`` | $(γmin) | ``1/s`` |
	| ``\dot{\gamma}_{\rm max}`` | $(γmax) | ``1/s`` |
	| ``N_{\dot{\gamma}}`` | $(Nγ) |  |
	| ``H`` | $(H) | ``m`` |
	| ``β`` | $(β) | ``Pa\cdot s / m`` |
	"""
end

# ╔═╡ e024bd4c-e501-42c8-9396-6900b3f5c583
md"""## Plotting"""

# ╔═╡ 35b29653-d0bb-46ee-9d28-c03ebf6a09e0
let
	τ = first.(solve_τ.(γrng; τ₁=τ₁, K=K, n=n, η₀=η₀, β=β, H=H))
	layout = plot(layout=(2,1))
	plot!(layout, γrng, τ; lw=2, color=1, label="", subplot=1, ylim=(0,maximum(τ)), ylabel=L"\tau~[Pa]")
	plot!(layout, γrng, τ; lw=2, color=1, label="", subplot=2, xscale=:log10, yscale=:log10, xlabel=L"\dot{\gamma}~[1/s]", ylabel=L"\tau~[Pa]")
end

# ╔═╡ db4ddd01-57a4-4746-9e0d-6ca3dff90fad
let
	if allow_slip
		layout = plot(layout=(2,1))
		for (i,logβ) in enumerate(5:8)
			τ = first.(solve_τ.(γrng; τ₁=τ₁, K=K, n=n, η₀=η₀, β=10.0^logβ, H=H))
			plot!(layout, γrng, τ; lw=2, color=i, label="10^$(logβ)", subplot=1, ylim=(0,maximum(τ)), legend=:bottomright, ylabel=L"\tau~[Pa]")
			plot!(layout, γrng, τ; lw=2, color=i, label="10^$(logβ)", subplot=2, xscale=:log10, yscale=:log10, legend=:bottomright, xlabel=L"\dot{\gamma}~[1/s]", ylabel=L"\tau~[Pa]")
		end
		plot!()
	end
end

# ╔═╡ Cell order:
# ╟─316a838d-57db-4194-b4e1-857f3afc44f6
# ╟─791fbd6c-1a96-4919-ac04-6ee5abc1a538
# ╟─2d4046cd-fd85-4f02-8126-469d296fda1a
# ╟─8c337acd-42f9-4850-8f6f-d3397a3183c1
# ╟─75987388-9b29-4505-b641-e1d4c8d72918
# ╟─dd4aebde-8087-402b-b17f-59b17522cda8
# ╟─99a62371-6939-4c1f-af94-eeb7d177ceb7
# ╟─f3407c7b-4686-4a9d-b354-5d3e9c8a4033
# ╟─8c15239f-c3e8-4e79-99a8-9b8c5b56bc25
# ╟─c58882b9-ca14-479f-ad4d-a98e194b5774
# ╟─a8ba3e65-e0a2-4a9f-bda4-5c4a4a64f7a1
# ╟─8c36449a-99f8-4bbc-9194-c432aae6708b
# ╟─2d129381-d823-403f-8eec-4322da51ddd7
# ╟─256361ad-327e-448c-ba72-cafccd4cea77
# ╟─2bf18333-273f-493c-93ae-0e0a9ae4eefe
# ╟─e37fa6c6-0fb5-41ec-9818-632f245131f7
# ╟─3772106b-1e2f-41a6-94bd-64a9459239c2
# ╟─3a3ead59-460a-4730-a012-c1244d163446
# ╟─0fc3ad2f-6742-4b08-9114-b8d831498463
# ╟─e024bd4c-e501-42c8-9396-6900b3f5c583
# ╟─35b29653-d0bb-46ee-9d28-c03ebf6a09e0
# ╟─db4ddd01-57a4-4746-9e0d-6ca3dff90fad
