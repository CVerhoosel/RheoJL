function solve_∂p∂r(h, τ₁, K, n, η₀, Q; ∂p∂r₀=0.0, rtol=1e-6, atol=1e-12, maxiter=100, bracket = [NaN, NaN], verbose=false)

    @assert Q isa Real "Q must be a real number"

    residual = ∂p∂r -> ( flux(h, ∂p∂r, τ₁, K, n, η₀) .- [Q,0.0] )

    if any(isnan.(bracket))
        ∂p∂r = -(3.0*η₀*Q)/(2.0*h^3)
        bracket[:] = sort([0., ∂p∂r])
    end

    ∂p∂r, info = newton(residual, x₀=∂p∂r₀, tol=rtol*abs(Q)+atol, maxiter=maxiter, verbose=verbose, bracket=bracket)

    if verbose
        for infoᵢ in info
            infoᵢ[2] += Q
        end
    end

	return ∂p∂r, info
end

function newton(residual; x₀=0.0, tol=1e-9, maxiter=25, bracket=[NaN, NaN], verbose=false)

    @assert x₀ isa Real "x₀ must be a real number"

    # Check the bracket is specified
    extrema = [NaN, NaN]
    if all(.!isnan.(bracket))
        @assert bracket isa AbstractVector{<:Real} && length(bracket) == 2 && bracket[1] < bracket[2] "Bracket must be a vector of two increasing real numbers."

        # Compute the extrema
        extrema[:] = first.(residual.(bracket))

        # Check whether an endpoint is a root
        for i in 1:2
            if abs(extrema[i]) < tol
                return bracket[i], verbose ? [[bracket[i], extrema[i], NaN, bracket[1], bracket[2]]] : 0
            end
        end

        @assert prod(extrema) < 0 "The bracket extrema should be of opposite sign."
        @assert bracket[1] <= x₀ <= bracket[2] "The initial guess x₀ should be within the bracket."
    end

    # Initialize the solution
    xᵢ = x₀
    rᵢ, ∂rᵢ = residual(xᵢ)

    # Convergence information
    if verbose
        info = [[xᵢ, rᵢ, ∂rᵢ, bracket[1], bracket[2]]]
    else
        info = []
    end

    # Check whether the initial solution has converged
    if abs(rᵢ) < tol
        return xᵢ, verbose ? info : 0
    end

    # Newton iterations
    for i in 1:maxiter

        # Store the previous iteration
        xᵢ₋₁ = xᵢ
        rᵢ₋₁ = rᵢ

        # Propose Newton update
        xₙ = xᵢ - rᵢ/∂rᵢ
        rₙ, ∂rₙ = residual(xₙ)

        if any(isnan.(bracket))
            # Accept the Newton proposal
            xᵢ = xₙ
            rᵢ = rₙ
            ∂rᵢ = ∂rₙ

            # Initialize the bracket
            if rᵢ₋₁*rᵢ < 0.0
                if xᵢ > xᵢ₋₁
                    bracket[:] = [xᵢ₋₁, xᵢ]
                    extrema[:] = [rᵢ₋₁, rᵢ]
                else
                    bracket[:] = [xᵢ, xᵢ₋₁]
                    extrema[:] = [rᵢ, rᵢ₋₁]
                end
            end
        else
            if bracket[1] < xₙ < bracket[2]
                # Accept the Newton proposal
                xᵢ = xₙ
                rᵢ = rₙ
                ∂rᵢ = ∂rₙ

                # Update the bracket
                if rᵢ*extrema[1] > 0.0
                    bracket[1] = xᵢ
                    extrema[1] = rᵢ
                else
                    bracket[2] = xᵢ
                    extrema[2] = rᵢ
                end
            else
                # Use bisection instead of the Newton proposal
                xᵢ = (bracket[1] + bracket[2]) / 2
                rᵢ, ∂rᵢ = residual(xᵢ)
            end
        end

        # Store the new iteration for output
        if verbose
            push!(info, [xᵢ, rᵢ, ∂rᵢ, bracket[1], bracket[2]])
        end

        # Check whether the solution has converged
        if abs(rᵢ) < tol
            if !verbose
                push!(info, i)
            end
            break
        end

        if i==maxiter
            error("Newton iteration did not converge in $(i) iterations")
        end
    end

	return xᵢ, verbose ? info : info[1]

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