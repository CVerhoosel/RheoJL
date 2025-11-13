function solve_τ(∂γ∂tₑ; τ₁=nothing, K=nothing, n=nothing, η₀=nothing, β=nothing, H=nothing, maxiter=20, tol=1e-10)
	@assert ∂γ∂tₑ ≥ 0. "Engineering γ̇ₑ shear rate should be non-negative"

	∂γ∂t = nothing
	if β===nothing && H===nothing
		∂γ∂t = ∂γ∂tₑ
	elseif β!==nothing && H!==nothing
		@assert β ≥ 0 "Slip parameter β must be non-negative"
		@assert H > 0 "Height H must be positive"
		∂γ∂t, τᵤ = solve_∂γ∂t(∂γ∂tₑ; τ₁=τ₁, K=K, n=n, η₀=η₀, βH=β*H, maxiter=maxiter, tol=tol)
	else
		throw("Both β and H must be specified if slip is incorporated.")
	end

	if η₀===nothing # Newtonian (K-variant), Bingham and Herschel-Bulkley
		@assert K > 0 "K must be positive"
		if n===nothing # Newtonian (K-variant) and Bingham
			if τ₁===nothing # Newtonian (K-variant)
				return K*∂γ∂t
			else # Bingham
				@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
				return τᵤ===nothing ? τ₁ + K*∂γ∂t : τᵤ
			end
		else # Herschel-Bulkley
			@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
			@assert n > 0 "n must be positive"
			return τᵤ === nothing ? τ₁ + K*∂γ∂t^n : τᵤ
		end
	else # Newtonian (η₀-variant) and biviscous (power law)
		@assert η₀ > 0 "Pre-yield viscosity η₀ must be positive"
		if n===nothing # Newtonian (η₀-variant) and biviscous
			if τ₁===nothing && K===nothing # Newtonian (η₀-variant)
				return η₀*∂γ∂t
			elseif τ₁!==nothing && K!==nothing #Biviscous
				@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
				@assert K > 0 "K must be positive"
				∂γ∂t₁ = τ₁/η₀
				τ₀ = τ₁ - K*∂γ∂t₁
				return ∂γ∂t < ∂γ∂t₁ ? η₀*∂γ∂t : τ₀ + K*∂γ∂t
			else
				if K===nothing
			  		throw("If η₀ and τ₁ are specified, so should be K.")
				else
					throw("If η₀ and K are specified, so should be τ₁.")
				end
			end
		else # Biviscous power law
			@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
			@assert K > 0 "K must be positive"
			@assert n > 0 "n must be positive"
			∂γ∂t₁ = τ₁/η₀
			τ₀ = τ₁ - K*∂γ∂t₁^n
			return ∂γ∂t < ∂γ∂t₁ ? η₀*∂γ∂t : τ₀ + K*∂γ∂t^n
		end
	end
end

function solve_∂γ∂t(∂γ∂tₑ; τ₁, K, n, η₀, βH, maxiter, tol)
	∂γ∂t = nothing	
	if η₀===nothing # Newtonian (K-variant), Bingham and Herschel-Bulkley
		@assert K > 0 "K must be positive"
		σ = K/βH # Dimensionless slip parameter
		if n===nothing # Newtonian (K-variant) and Bingham
			if τ₁===nothing # Newtonian (K-variant)
				return ∂γ∂tₑ/(1 + σ), nothing
			else # Bingham
				@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
				return βH*∂γ∂tₑ < τ₁ ? (0., βH*∂γ∂tₑ) : ((∂γ∂tₑ - τ₁/βH)/(1 + σ), nothing)
			end
		else # Herschel-Bulkley
			@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
			@assert n > 0 "n must be positive"
			if βH*∂γ∂tₑ < τ₁
				return 0., βH*∂γ∂tₑ
			else # Yielding
				r = x -> (x - ∂γ∂tₑ + τ₁/βH + σ*x^n, 1 + n*σ*x^(n-1), [])
				∂γ∂t, info = newton(r, x₀=tol, bracket=[0., Inf], maxiter=maxiter, tol=tol)
				return ∂γ∂t, nothing
			end
		end
	else # Newtonian (η₀-variant) and biviscous (power law)
		@assert η₀ > 0 "Pre-yield viscosity η₀ must be positive"
		μ = η₀/βH # Dimensionless slip parameter
		if n===nothing # Newtonian (η₀-variant) and biviscous
			if τ₁===nothing && K===nothing # Newtonian (η₀-variant)
				return ∂γ∂tₑ/(1 + μ), nothing
			elseif τ₁!==nothing && K!==nothing #Biviscous
				@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
				@assert K > 0 "K must be positive"
				∂γ∂t₁ = τ₁/η₀
				τ₀ = τ₁ - K*∂γ∂t₁
				return (μ/(1 + μ))*βH*∂γ∂tₑ < τ₁ ? ∂γ∂tₑ/(1 + μ) : (∂γ∂tₑ - τ₀/βH)/(1 + K/βH), nothing
			else
				if K===nothing
			  		throw("If η₀ and τ₁ are specified, so should be K.")
				else
					throw("If η₀ and K are specified, so should be τ₁.")
				end
			end
		else # Biviscous power law
			@assert τ₁ ≥ 0 "Yield stress τ₁ must be non-negative"
			@assert K > 0 "K must be positive"
			@assert n > 0 "n must be positive"
			∂γ∂t₁ = τ₁/η₀
			τ₀ = τ₁ - K*∂γ∂t₁^n
			if (μ/(1 + μ))*βH*∂γ∂tₑ < τ₁
				return ∂γ∂tₑ/(1 + μ), nothing
			else # Yielding				
				r = x -> (x - ∂γ∂tₑ + τ₀/βH + (K/βH)*x^n, 1 + n*(K/βH)*x^(n-1), [])
				∂γ∂t, info = newton(r, x₀=∂γ∂t₁, bracket=[∂γ∂t₁, Inf], maxiter=maxiter, tol=tol)
				return ∂γ∂t, nothing
			end
		end
	end
end