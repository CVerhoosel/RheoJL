function residual_∂p∂r(h, ∂p∂r, τ₁, K, n, η₀, Q; β=nothing)

    Qᵢ, ∂Q∂∂p∂rᵢ = flux(h, ∂p∂r, τ₁, K, n, η₀; β=β)
    
    r       = Qᵢ - Q
    ∂r∂∂p∂r = ∂Q∂∂p∂rᵢ

    return r, ∂r∂∂p∂r, []
end

function solve_∂p∂r(h, τ₁, K, n, η₀, Q; β=nothing, ∂p∂r₀=0.0, rtol=1e-6, atol=1e-12, maxiter=1_000, bracket=[-Inf, Inf], output=:short)

    @assert Q isa Real && Q ≥ -eps()  "Q must be a non-negative real number, since this should be enforced via the outer loop bracket."

    residual = ∂p∂r -> residual_∂p∂r(h, ∂p∂r, τ₁, K, n, η₀, Q; β=β)

    ∂p∂r, info = newton(residual; x₀=∂p∂r₀, tol=rtol*abs(Q)+atol, maxiter=maxiter, bracket=bracket, output=output)

    if output==:long
        for infoᵢ in info
            infoᵢ[2] += Q
        end
    end

	return ∂p∂r, info
end

function residual_force(∂h∂t, h, τ₁, K, n, η₀, F, rᵥ; β=nothing, γ=nothing, α=nothing, ∂p∂r₀ₘ=nothing, output=:short, rtolinner=1e-6, atolinner=1e-12, maxiterinner=1_000, bracketinner=[-Inf, Inf], outputinner=:short)

    @assert ∂h∂t isa Real && ∂h∂t ≤ eps() "∂h∂t must be a non-positive real number, since this should be enforced via the outer loop bracket."
    @assert all(rᵥ .≥ -eps()) && all(isreal, rᵥ) "All rᵥ must be non-negative real numbers"
    
    rₘ = vertex_to_midpoint(rᵥ)
    R  = rᵥ[end]

    # Initialize pressure gradient estimate
    if ∂p∂r₀ₘ === nothing
        ∂p∂r₀ₘ = -(4*F)/(π*R^4) * rₘ
    end

    # Evaluate the flux
    Qₘ = -∂h∂t*rₘ

	results = [solve_∂p∂r(h, τ₁, K, n, η₀, Qₘ[i]; β=β, ∂p∂r₀=∂p∂r₀ₘ[i], rtol=rtolinner, atol=atolinner, maxiter=maxiterinner, bracket=bracketinner, output=outputinner) for i in eachindex(Qₘ)]
    ∂p∂rₘ = getindex.(results, 1)
    info  = getindex.(results, 2)

    if output==:long
        if outputinner==:long
            iterations = size.(info,1) .- 1
        else
            iterations = info
        end
    else
        iterations = []
    end

    ∂Q∂∂p∂rₘ = last.(flux.(h, ∂p∂rₘ, τ₁, K, n, η₀; β=β))
    ∂∂p∂r∂∂h∂tₘ = -rₘ ./ ∂Q∂∂p∂rₘ

    # Compute the residual
    r       = force(rᵥ, ∂p∂rₘ; h=h, γ=γ, α=α) - F
    ∂r∂∂h∂t = force(rᵥ, ∂∂p∂r∂∂h∂tₘ; h=h, γ=γ, α=α)

    return r, ∂r∂∂h∂t, iterations
end

function solve_system(h, R, τ₁, K, n, η₀, F, N; β=nothing, γ=nothing, α=nothing, ∂h∂t₀=0.0, ∂p∂r₀ₘ=nothing, rtol=1e-6, atol=1e-12, maxiter=25, output=:short, rtolinner=1e-6, atolinner=1e-12, maxiterinner=1_000, bracketinner=[-Inf, Inf], outputinner=:short)

    # Get the vertices mesh
    rᵥ = collect(range(0, R, length=N+1))

    # Formulate the residual
    residual = ∂h∂t -> residual_force(∂h∂t, h, τ₁, K, n, η₀, F, rᵥ; β=β, γ=γ, α=α, ∂p∂r₀ₘ=∂p∂r₀ₘ, output=output, rtolinner=rtolinner, atolinner=atolinner, maxiterinner=maxiterinner, bracketinner=bracketinner, outputinner=outputinner)

    # Prevent positive height rates
    bracket=[-Inf, 0]

    # Perform Newton iterations on the force residual
    ∂h∂t, info = newton(residual, x₀=∂h∂t₀, tol=rtol*abs(F)+atol, maxiter=maxiter, bracket=bracket, output=output)

    if output==:long
        for infoᵢ in info
            infoᵢ[2] += F
        end
    end

	return ∂h∂t, info
end

function newton(residual; x₀=0.0, tol=1e-9, maxiter=25, bracket=[-Inf, Inf], output=:short)

    @assert output in [:short, :long] "Output must be either :short or :long."
    @assert bracket isa AbstractVector{<:Real} && length(bracket) == 2 && bracket[1] < bracket[2] "Bracket must be a vector of two increasing real numbers."
    @assert x₀ isa Real && isfinite(x₀) && bracket[1] <= x₀ <= bracket[2] "x₀ must be a finite real number inside or on the boundary of the bracket"

    # Copy the bracket to avoid by-reference manipulation
    bracket = copy(bracket)

    # Evaluate the extrema
    extrema = [-Inf, Inf]
    for i in 1:2
        if isfinite(bracket[i])
            extrema[i], _, iterstats = residual(bracket[i])
        end

        if abs(extrema[i]) < tol
            return bracket[i], output==:long ? [[bracket[i], extrema[i], bracket[1], bracket[2], iterstats...]] : 0
        end
    end

    # Correct the unbounded interval extremum sign if required
    if isfinite(extrema[1]) && isinf(extrema[2]) && prod(extrema) > 0
        extrema[2] *= -1
    elseif isinf(extrema[1]) && isfinite(extrema[2]) && prod(extrema) > 0
        extrema[1] *= -1
    end

    @assert prod(extrema) < 0 "The bracket extrema should be of opposite sign."
    
    # Initialize the solution
    xᵢ = copy(x₀)
    rᵢ, ∂rᵢ, iterstats = residual(xᵢ)

    # Convergence information
    info = [[xᵢ, rᵢ, bracket[1], bracket[2], iterstats...]]

    # Check whether the initial solution has converged
    if abs(rᵢ) < tol
        return xᵢ, output==:long ? info : 0
    end

    # Newton iterations
    for i in 1:maxiter

        # Store the previous iteration
        xᵢ₋₁ = xᵢ
        rᵢ₋₁ = rᵢ

        # Propose Newton update
        xₙ = xᵢ - rᵢ/∂rᵢ
        
        if all(isinf.(bracket))
            # Bracket is unbounded on both sides, i.e., [-∞,∞]

            # Accept the Newton proposal
            xᵢ = xₙ
            rᵢ, ∂rᵢ, iterstats = residual(xᵢ)

            # Initialize the bracket
            if rᵢ*rᵢ₋₁ < 0.0
                if xᵢ > xᵢ₋₁
                    bracket = [xᵢ₋₁, xᵢ]
                    extrema = [rᵢ₋₁, rᵢ]
                else
                    bracket = [xᵢ, xᵢ₋₁]
                    extrema = [rᵢ, rᵢ₋₁]
                end
            end
        else
            if isinf(bracket[1])
                # Bracket is unbounded on the left, i.e., [-∞,right]
                if xₙ < bracket[2]
                    xᵢ = xₙ
                else
                    xᵢ = bracket[2]
                end
            elseif isinf(bracket[2])
                # Bracket is unbounded on the right, i.e., [left, ∞]
                if xₙ > bracket[1]
                    xᵢ = xₙ
                else
                    xᵢ = bracket[1]
                end
            else
                # Bracket is bounded on both sides, i.e., [left, right]
                if bracket[1] < xₙ < bracket[2]
                    xᵢ = xₙ
                else
                    # Use bisection instead of the Newton proposal if the proposal is not within the bracket
                    xᵢ = (bracket[1] + bracket[2]) / 2
                end
            end

            # Update the bracket
            rᵢ, ∂rᵢ, iterstats = residual(xᵢ)

            if rᵢ*extrema[1] > 0.0
                bracket[1] = xᵢ
                extrema[1] = rᵢ
            else
                bracket[2] = xᵢ
                extrema[2] = rᵢ
            end
        end

        # Store the new iteration for output
        if output==:long
            push!(info, [xᵢ, rᵢ, bracket[1], bracket[2], iterstats...])
        end

        # Check whether the solution has converged
        if abs(rᵢ) < tol
            return xᵢ, output==:long ? info : i
        end
    end

    error("Newton solver did not converge in $(maxiter) iterations")
end

function integrate_system(h₀, R₀, τ₁, K, n, η₀, F, Nᵣ, T, Nₜ, Δt₀; β=nothing, γ=nothing, α=nothing, Δₘₐₓ=0.1, target=5,
    rtol=1e-6, atol=1e-12, maxiter=100, output=:short, rtolinner=1e-6, atolinner=1e-12, maxiterinner=1_000, outputinner=:short)
 
    @assert output in [:short, :long] "Output must be either :short or :long."

    # Initialization
    t    = 0.0
	∂h∂t = 0.0
    R    = R₀
    h    = h₀
    V    = 2*h*π*R^2
    ∂R∂t = -R*∂h∂t/(2*h)

    # Define the initial pressure gradient
    rᵥ    = collect(range(0, R₀, length=Nᵣ+1))
    rₘ    = vertex_to_midpoint(rᵥ)
    ∂p∂rₘ = -(4*F)/(π*R^4) * rₘ

    # Define the geometric time series
    tₛ = geometric_time_sequence(T, Δt₀, Nₜ)
    Δtₛ = tₛ[2:end] - tₛ[1:end-1]
    ρ = Δtₛ[end] / Δtₛ[end-1]
    Δtₛ = vcat(Δtₛ, [ρ*Δtₛ[end]])
    Δtₛ = linear_interpolation(tₛ, Δtₛ)

    # Initialize the time step scale
    σ = 1.0

    # Initialize the solution and info storage
	sol  = [[t σ R]]
    info = []

    while T-t > atol

        # Update the height rate
        ∂h∂t, ∂h∂t_info = solve_system(h, R, τ₁, K, n, η₀, F, Nᵣ; β=β, γ=γ, α=α, ∂h∂t₀=∂h∂t, ∂p∂r₀ₘ=∂p∂rₘ,
            rtol=rtol, atol=atol, maxiter=maxiter, output=output, 
            rtolinner=rtolinner, atolinner=atolinner, maxiterinner=maxiterinner, outputinner=outputinner)

        info = push!(info, ∂h∂t_info)

        # Update the pressure gradient
        Qₘ      = -∂h∂t*rₘ
        results = [solve_∂p∂r(h, τ₁, K, n, η₀, Qₘ[i]; β=β, ∂p∂r₀=∂p∂rₘ[i], rtol=rtolinner, atol=atolinner, maxiter=maxiterinner, output=outputinner) for i in eachindex(Qₘ)]
        ∂p∂rₘ   = getindex.(results, 1)

        # Get the current maximum time step as per the geometric series
        Δtₘₐₓ = Δtₛ(t)

        # Adjust the time step scale based on the number of Newton iterations
        if target !== nothing
            σ = min(1.0, σ*(target/ (output==:short ? ∂h∂t_info+1 : size(∂h∂t_info,1))))
        end

        # Scale the time step
        Δt = σ*Δtₘₐₓ

        # Limit the time step scale based on the maximum allowed height change
        if Δₘₐₓ !== nothing && Δt*∂h∂t < -Δₘₐₓ*h
            Δt = -Δₘₐₓ*h/∂h∂t
        end
        
        # Clip the time step to not exceed the final time
        if Δt > T-t
            Δt = T-t
        end

        # Update the height, radius and time
        t    = t + Δt
        h    = h + ∂h∂t * Δt
        R    = sqrt(V/(2*π*h))
        ∂R∂t = -R*∂h∂t/(2*h)
        
        # Store the solution
        push!(sol, [t σ R])
    end

    return reduce(vcat, sol), info
end