function velocity_profile(h, ∂p∂r, τ₁, K, n, η₀, np=1_000)

    @assert τ₁ ≥ 0 "τ₁ must be non-negatve"
    @assert η₀ > 0 "η₀ must be positive" 
    @assert n > 0 "n must be positive"
    @assert K > 0 "K must be positive"
 
    z = collect(0:np-1) / (np-1) * h

    if τ₁ ≥ h * abs(∂p∂r)
        v = 0.5*∂p∂r/η₀ * (z.^2 .- h^2)
        np1 = np
    else
        wy = τ₁ / abs(∂p∂r)

        τ₀ = τ₁ - K * (τ₁/η₀)^n

        np1 = sum(z .< wy)
        
        α = ∂p∂r / K
        β = τ₀ / K
        ρ = (n+1)/n

        # Centerline velocity
        vc = (1/(α*ρ)) * ((abs(α)*wy - β)^ρ - (abs(α)*h - β)^ρ) - 0.5 * (∂p∂r / η₀) * wy^2

        v = 0.5*∂p∂r/η₀ * z.^2 .+ vc

        z2 = z[np1+1:end]
        v2 = (1/(α*ρ)) * ((abs(α)*z2 .- β).^ρ .- (abs(α)*h .- β).^ρ)

        v[np1+1:end] = v2
    end

    return z, v, np1
end

function flux(h, ∂p∂r, τ₁, K, n, η₀)

    @assert τ₁ ≥ 0 "τ₁ must be non-negatve"
    @assert η₀ > 0 "η₀ must be positive" 
    @assert n  > 0 "n must be positive"
    @assert K  > 0 "K must be positive"
    @assert h  > 0 "h must be positive"
 

    if τ₁ ≥ h * abs(∂p∂r)
        Q1 = -2.0*(∂p∂r*h^3)/(3.0*η₀)
        Q2 = 0.0

        ∂Q1 = -2.0*h^3/(3.0*η₀)
        ∂Q2 = 0.0
    else
        τ₀ = τ₁ - K * (τ₁/η₀)^n

        wy = τ₁ / abs(∂p∂r)
        ∂wy = -sign(∂p∂r) * τ₁ / (abs(∂p∂r)^2) 
        
        α  = ∂p∂r / K
        ∂α = 1 / K
        β  = τ₀ / K
        ρ  = (n+1)/n

        # Centerline velocity
        vc  = (1/(α*ρ)) * ((abs(α)*wy - β)^ρ - (abs(α)*h - β)^ρ) - 0.5 * (∂p∂r / η₀) * wy^2
        ∂vc = -(1/(α^2*ρ)) * ∂α * ((abs(α)*wy - β)^ρ - (abs(α)*h - β)^ρ) + 
              (1/α) * ((sign(α)*∂α*wy+abs(α)*∂wy)*(abs(α)*wy - β)^(ρ-1) - (sign(α)*∂α*h)*(abs(α)*h - β)^(ρ-1)) +
              - 0.5 * (wy^2) / η₀ +
              - (∂p∂r / η₀) * wy * ∂wy

        Q1 = (∂p∂r*wy^3)/(3.0*η₀) + 2.0*wy*vc
        ∂Q1 = (wy^3)/(3.0*η₀) + (3*∂p∂r*wy^2*∂wy)/(3.0*η₀) + 2.0*∂wy*vc + 2.0*wy*∂vc

        Q2  = (2.0/(α*abs(α)*ρ*(ρ+1))) * ((abs(α)*h - β)^(ρ+1) - (abs(α)*wy - β)^(ρ+1)) - (2.0*(h-wy)/(ρ*α))*(abs(α)*h - β)^ρ
        ∂Q2 = -((4.0*∂α)/((α^2)*abs(α)*ρ*(ρ+1))) * ((abs(α)*h - β)^(ρ+1) - (abs(α)*wy - β)^(ρ+1)) +
              (2.0/(α*abs(α)*ρ)) * ((sign(α)*∂α*h)*(abs(α)*h - β)^(ρ) - (sign(α)*∂α*wy+abs(α)*∂wy)*(abs(α)*wy - β)^(ρ)) +
              (2.0*(∂wy)/(ρ*α))*(abs(α)*h - β)^ρ +
              (2.0*(h-wy)/(ρ*α^2))*∂α*(abs(α)*h - β)^ρ +
              - (2.0*(h-wy)/(α))*(sign(α)*∂α*h)*(abs(α)*h - β)^(ρ-1)
    end

    return Q1+Q2, ∂Q1+∂Q2 
end

function force(rᵥ, ∂p∂rₘ)
	p = right_integrate(rᵥ, ∂p∂rₘ)
	F = integral(rᵥ, vertex_to_midpoint(2*π*rᵥ.*p))
	return F
end