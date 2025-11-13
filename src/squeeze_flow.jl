"""
    shear_stress(γ, τ₁, K, n, η₀)

Compute the shear stress using the regularized Herschel-Bulkley model.

The regularized Herschel-Bulkley model provides a smooth transition from 
Newtonian behavior at low shear rates to power-law behavior at high shear rates,
avoiding the discontinuity at the yield stress.

# Arguments
- `γ`: Shear rate (1/s)
- `τ₁`: Yield stress (Pa)
- `K`: Consistency index (Pa⋅sⁿ)
- `n`: Flow behavior index (dimensionless)
- `η₀`: Newtonian viscosity at low shear rates (Pa⋅s)

# Returns
- `τ`: Shear stress (Pa)

# Model description
For low shear rates (η₀|γ| > τ₁):
    τ = η₀γ  (Newtonian behavior)

For high shear rates (η₀|γ| ≤ τ₁):
    τ = sign(γ)(τ₀ + K|γ|ⁿ)  (Power-law behavior)
    
where τ₀ = τ₁ - K(τ₁/η₀)ⁿ ensures continuity at the transition point.
"""
function shear_stress(γ, τ₁, K, n, η₀)
    @assert τ₁ ≥ 0 "τ₁ must be non-negatve"
    @assert η₀ > 0 "η₀ must be positive" 
    @assert n > 0 "n must be positive"
    @assert K > 0 "K must be positive"

    if η₀ * abs(γ) ≤ τ₁
        τ = η₀*γ
    else
        τ₀ = τ₁ - K * (τ₁/η₀)^n
        τ = sign(γ)*(τ₀ + K*abs(γ)^n)
    end

    return τ
end

function velocity_profile(h, ∂p∂r, τ₁, K, n, η₀; β=nothing, np=1_000)

    @assert τ₁ ≥ 0 "τ₁ must be non-negatve"
    @assert η₀ > 0 "η₀ must be positive" 
    @assert n > 0 "n must be positive"
    @assert K > 0 "K must be positive"
    
    # Slip velocity
    vₛ = 0.0
    if β !== nothing
        @assert β > 0 "β must be positive"
        τₕ = ∂p∂r * h
        vₛ = - τₕ / β
    end
 
    z  = collect(0:np-1) / (np-1) * h

    if τ₁ ≥ h * abs(∂p∂r)
        v = 0.5*∂p∂r/η₀ * (z.^2 .- h^2) .+ vₛ
        np1 = np
    else
        wy = τ₁ / abs(∂p∂r)

        τ₀ = τ₁ - K * (τ₁/η₀)^n

        np1 = sum(z .< wy)
        
        α = ∂p∂r / K
        γ = τ₀ / K
        ρ = (n+1)/n

        # Centerline velocity
        vc = (1/(α*ρ)) * ((abs(α)*wy - γ)^ρ - (abs(α)*h - γ)^ρ) - 0.5 * (∂p∂r / η₀) * wy^2

        v = 0.5* (∂p∂r / η₀) * z.^2 .+ vc .+ vₛ

        z2 = z[np1+1:end]
        v2 = (1/(α*ρ)) * ((abs(α)*z2 .- γ).^ρ .- (abs(α)*h .- γ).^ρ) .+ vₛ

        v[np1+1:end] = v2
    end

    return z, v, np1
end

function flux(h, ∂p∂r, τ₁, K, n, η₀; β=nothing)

    @assert τ₁ ≥ 0 "τ₁ must be non-negatve"
    @assert η₀ > 0 "η₀ must be positive" 
    @assert n  > 0 "n must be positive"
    @assert K  > 0 "K must be positive"
    @assert h  > 0 "h must be positive"

    # Slip velocity and slip
    vₛ = 0.0
    Q₃ = 0.0
    ∂Q₃ = 0.0
    if β !== nothing
        @assert β > 0 "β must be positive"
        τₕ = ∂p∂r * h
        vₛ = - τₕ / β
        Q₃ = 2.0*h*vₛ
        ∂Q₃ = -2.0*h^2/β
    end

    if τ₁ ≥ h * abs(∂p∂r)
        Q₁ = -2.0*(∂p∂r*h^3)/(3.0*η₀)
        Q₂ = 0.0
        
        ∂Q₁ = -2.0*h^3/(3.0*η₀)
        ∂Q₂ = 0.0
    else
        τ₀ = τ₁ - K * (τ₁/η₀)^n

        wy = τ₁ / abs(∂p∂r)
        ∂wy = -sign(∂p∂r) * τ₁ / (abs(∂p∂r)^2) 
        
        α  = ∂p∂r / K
        ∂α = 1 / K
        γ  = τ₀ / K
        ρ  = (n+1)/n

        # Centerline velocity
        vc  = (1/(α*ρ)) * ((abs(α)*wy - γ)^ρ - (abs(α)*h - γ)^ρ) - 0.5 * (∂p∂r / η₀) * wy^2
        ∂vc = -(1/(α^2*ρ)) * ∂α * ((abs(α)*wy - γ)^ρ - (abs(α)*h - γ)^ρ) + 
              (1/α) * ((sign(α)*∂α*wy+abs(α)*∂wy)*(abs(α)*wy - γ)^(ρ-1) - (sign(α)*∂α*h)*(abs(α)*h - γ)^(ρ-1)) +
              - 0.5 * (wy^2) / η₀ +
              - (∂p∂r / η₀) * wy * ∂wy

        Q₁ = (∂p∂r*wy^3)/(3.0*η₀) + 2.0*wy*vc
        ∂Q₁ = (wy^3)/(3.0*η₀) + (3*∂p∂r*wy^2*∂wy)/(3.0*η₀) + 2.0*∂wy*vc + 2.0*wy*∂vc

        Q₂  = (2.0/(α*abs(α)*ρ*(ρ+1))) * ((abs(α)*h - γ)^(ρ+1) - (abs(α)*wy - γ)^(ρ+1)) - (2.0*(h-wy)/(ρ*α))*(abs(α)*h - γ)^ρ
        ∂Q₂ = -((4.0*∂α)/((α^2)*abs(α)*ρ*(ρ+1))) * ((abs(α)*h - γ)^(ρ+1) - (abs(α)*wy - γ)^(ρ+1)) +
              (2.0/(α*abs(α)*ρ)) * ((sign(α)*∂α*h)*(abs(α)*h - γ)^(ρ) - (sign(α)*∂α*wy+abs(α)*∂wy)*(abs(α)*wy - γ)^(ρ)) +
              (2.0*(∂wy)/(ρ*α))*(abs(α)*h - γ)^ρ +
              (2.0*(h-wy)/(ρ*α^2))*∂α*(abs(α)*h - γ)^ρ +
              - (2.0*(h-wy)/(α))*(sign(α)*∂α*h)*(abs(α)*h - γ)^(ρ-1)
    end

    return Q₁+Q₂+Q₃, ∂Q₁+∂Q₂+∂Q₃
end

function force(rᵥ, ∂p∂rₘ; h=nothing, γ=nothing, α=nothing)

    # Compute the laplace pressure
    Δp = 0.0
    if γ !== nothing && α !== nothing
        @assert γ isa Real && γ ≥ 0 "γ must be a non-negative real number"
        @assert α isa Real && 0 ≤ α ≤ 1 "α must be between 0 and 1"
        @assert h > 0 "h must be positive to compute the Laplace pressure"
        Δp = 2*γ*α/h
    elseif γ !== nothing || α !== nothing
        error("Both γ (surface tension) and α (relative curvature) must be specified for capillary effects.")
    end

	p = right_integrate(rᵥ, ∂p∂rₘ) .+ Δp
	F = integral(rᵥ, vertex_to_midpoint(2*π*rᵥ.*p))

	return F
end