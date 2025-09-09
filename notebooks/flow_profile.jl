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
	
	using Revise, PlutoUI, Plots, LaTeXStrings
	using RheoJL

md"""
!!! info "Julia environment"
	This notebook makes use of the global environment. Make sure all used packages, including RheoJL, are available there.
"""
end

# ╔═╡ 316a838d-57db-4194-b4e1-857f3afc44f6
md"""
# RheoJL flow profile test
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
	``\tau_{y}~[Pa]``: $(@bind τ₁ Slider(0:1:200; default=130, show_value=true))
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

# ╔═╡ 3a3ead59-460a-4730-a012-c1244d163446
begin
	reset_parameters
	md"``R_0~[mm]``: $(@bind R₀_input Slider(5:0.1:25; default=13, show_value=true))"
end

# ╔═╡ 584713b9-455a-4df4-9694-de8acf84f801
begin
	reset_parameters
	md"""
	``V~[ml]``: $(@bind V_input Slider(0.1:0.01:1.5; default=0.73, show_value=true))
	"""
end

# ╔═╡ 2d129381-d823-403f-8eec-4322da51ddd7
begin
	V  = V_input / 1e6
	R₀ = R₀_input / 1e3
	h  = (V/(π*R₀^2))/2
	md"""
	| Parameter | Value | Unit |
	|-----------|-------|------|
	| ``R_0`` | $(R₀) | ``m^3`` |
	| ``V`` | $(V) | ``m^3`` |
	| ``h`` | $(h) | ``m`` |
	"""
end

# ╔═╡ e024bd4c-e501-42c8-9396-6900b3f5c583
md"""## Plotting
Set the pressure gradient:
"""

# ╔═╡ 2088b8d6-d9b1-4f8d-bcd8-44f872b333a9
let
	reset_parameters
	dpdr_max = 4*τ₁/h
	dpdr_range = range(-dpdr_max,dpdr_max,length=100)
	md"""
	``\frac{\partial p}{\partial r}``: $(@bind dpdr Slider(dpdr_range; default=0.5*dpdr_max, show_value=true))
	"""
end

# ╔═╡ f68edc63-fdf9-497a-89cb-355ca91032e6
z, v, nw = velocity_profile(h, dpdr, τ₁, K, n, η₀);

# ╔═╡ d5788ee3-73a0-432d-b7f3-eeb6f7c5b706
let
	∂γ∂t_max = abs((v[end]-v[end-1])/(z[end]-z[end-1]))
	plot([-∂γ∂t₀,∂γ∂t₀], [-τ₁,τ₁], linewidth=3, xlabel=L"\dot{\gamma}~[1/s]", ylabel=L"\tau~[Pa]", label="No yielding")
	if ∂γ∂t_max > ∂γ∂t₀
		∂γ∂t = range(∂γ∂t₀, ∂γ∂t_max, length=100)
		plot!(∂γ∂t, τ₀.+K*∂γ∂t.^n; linewidth=3, color=2, label="Yielding")
		plot!(-∂γ∂t, -τ₀.-K*∂γ∂t.^n; linewidth=3, color=2, label="")
	end
	plot!(title="Constitutive behavior")
end

# ╔═╡ 44d74c6c-5618-4290-84e3-a48a957a752c
let
	plot(v[1:nw], z[1:nw], label="Non-yielding", linewidth=3, xlabel=L"v~[m/s]", ylabel=L"z~[m]")
	plot!(v[nw+1:end], z[nw+1:end], label="Yielding", linewidth=3, title="Velocity profile")
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
# ╟─3a3ead59-460a-4730-a012-c1244d163446
# ╟─584713b9-455a-4df4-9694-de8acf84f801
# ╟─e024bd4c-e501-42c8-9396-6900b3f5c583
# ╟─2088b8d6-d9b1-4f8d-bcd8-44f872b333a9
# ╟─f68edc63-fdf9-497a-89cb-355ca91032e6
# ╟─d5788ee3-73a0-432d-b7f3-eeb6f7c5b706
# ╟─44d74c6c-5618-4290-84e3-a48a957a752c
