function solve_∂p∂r(h, τ₁, K, n, η₀, Q; ∂p∂r₀=0.0, rtol=1e-6, atol=1e-12, maxiter=100, verbose=false)

    @assert Q isa Real "Q must be a real number"

    residual = ∂p∂r -> ( flux(h, ∂p∂r, τ₁, K, n, η₀) .- [Q,0.0] )

    ∂p∂r, info = newton(residual, x₀=∂p∂r₀, tol=rtol*abs(Q)+atol, maxiter=maxiter, verbose=verbose)

    if verbose
        for infoᵢ in info
            infoᵢ[2] += Q
        end
    end

	return ∂p∂r, info
end

function newton(residual; x₀=0.0, tol=1e-9, maxiter=25, verbose=false)

    @assert x₀ isa Real "x₀ must be a real number"

    # Convergence information
    info = []

    # Initialize the solution
    x = x₀

    # Newton iterations
    for i in 1:maxiter

        # Compute the flux and its derivative w.r.t. ∂p∂r
        rᵢ, ∂rᵢ = residual(x)

        if verbose
            push!(info, [x, rᵢ, ∂rᵢ])
        end

        # Check whether the solution has converged
        if abs(rᵢ) < tol
            if !verbose
                push!(info, i)
            end

            break
        else
            # Update the solution
            x -= rᵢ/∂rᵢ
        end

        if i==maxiter
            error("Newton iteration did not converge in $(i) iterations")
        end
    end

	return x, verbose ? info : info[1]

end

# function eval_force(rᵥ, ∂p∂rₘ)
# 	p = right_integrate(rᵥ, ∂p∂rₘ)
# 	F = integral(rᵥ, vertex_to_midpoint(2*π*rᵥ.*p))
# end

# function solve_newton(h, R₀, τ₁, K, n, η₀, F, N; dhdt₀=0.0, rtol=1e-6, maxouter=100)
# 	rᵥ = collect(range(0, R₀, length=N))
# 	rₘ = vertex_to_midpoint(rᵥ)
# 	allinfo = zeros(length(rₘ),0)
	
# 	# Initialize dhdt and ∂p∂rₘ
# 	dhdt  = dhdt₀
# 	∂p∂rₙ = -(4*F)/(π*R₀^4)*rₘ 

# 	# Newton iterations
# 	for iiter in 1:maxouter

# 		# Solve for ∂p∂rₘ given dhdt
# 		∂p∂rₙ, info = solve_∂p∂rₘ(h, τ₁, K, n, η₀, rₘ, ∂p∂rₙ, dhdt)
# 		allinfo = hcat(allinfo, info)
		
# 		# Update dhdt
# 		∂Q = last.(Flux.(h, ∂p∂rₙ, τ₁, K, n, η₀))
	
# 		∂²p∂r∂hdotₙ  = -rₘ ./ ∂Q
# 		∂F∂dhdt = eval_force(rᵥ, ∂²p∂r∂hdotₙ)

# 		Fₙ   = eval_force(rᵥ, ∂p∂rₙ)
# 		dhdt = dhdt + (F-Fₙ)/∂F∂dhdt

# 		if abs(F-Fₙ) < rtol*F
# 			println("Newton algorithm converged in $(iiter) iterations.")
# 			return dhdt, Fₙ, ∂p∂rₙ, allinfo
# 		end
# 	end

# 	error("Newton algorithm did not converge in $(maxouter) iterations")
# end

# function integrate(h, R, τ₁, K, n, η₀, F, N, T, Δt₀; targetiter=6)
# 	t = 0.0
# 	Δt = Δt₀
# 	dhdt = 0.0
# 	V = 2*h*π*R^2

# 	sol = [t h R]
# 	allinfo = [NaN]
# 	while t < T
# 		if t + Δt ≥ T
# 			Δt = T - sol[end,1]
# 		end
		
# 		dhdt, F, ∂p∂r, info = solve_newton(h, R, τ₁, K, n, η₀, F, N, dhdt₀=dhdt)

# 		push!(allinfo, size(info,2))
		
# 		h = h + dhdt * Δt
# 		R = sqrt(V/(2*π*h))
# 		t = t + Δt

# 		sol = vcat(sol, [t h R])
		
# 		Δt = Δt * targetiter / size(info,2)
# 	end
	
# 	return sol, allinfo
# end