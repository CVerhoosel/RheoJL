function velocity_profile(h, dpdz, τy, K, n, η0, np=1_000)

    @assert τy ≥ 0 "τy must be non-negatve"
    @assert η0 > 0 "η0 must be positive" 
    @assert n > 0 "n must be positive"
    @assert K > 0 "K must be positive"
 
    z = collect(0:np-1) / (np-1) * h

    if τy ≥ h * abs(dpdz)
        v = 0.5*dpdz/η0 * (z.^2 .- h^2)
        np1 = np
    else
        wy = τy / abs(dpdz)

        τ0 = τy - K * (τy/η0)^n

        np1 = sum(z .< wy)
        
        α = dpdz / K
        β = τ0 / K
        ρ = (n+1)/n

        # Centerline velocity
        vc = (1/(α*ρ)) * ((abs(α)*wy - β)^ρ - (abs(α)*h - β)^ρ) - 0.5 * (dpdz / η0) * wy^2

        v = 0.5*dpdz/η0 * z.^2 .+ vc

        z2 = z[np1+1:end]
        v2 = (1/(α*ρ)) * ((abs(α)*z2 .- β).^ρ .- (abs(α)*h .- β).^ρ)

        v[np1+1:end] = v2
    end

    return z, v, np1
end

function flux(h, dpdz, τy, K, n, η0)

    @assert τy ≥ 0 "τy must be non-negatve"
    @assert η0 > 0 "η0 must be positive" 
    @assert n > 0 "n must be positive"
    @assert K > 0 "K must be positive"
    @assert h > 0 "h must be positive"
 

    if τy ≥ h * abs(dpdz)
        Q1 = -2.0*(dpdz*h^3)/(3.0*η0)
        Q2 = 0.0

        ∂Q1 = -2.0*h^3/(3.0*η0)
        ∂Q2 = 0.0
    else
        τ0 = τy - K * (τy/η0)^n

        wy = τy / abs(dpdz)
        ∂wy = -sign(dpdz) * τy / (abs(dpdz)^2) 
        
        α  = dpdz / K
        ∂α = 1 / K
        β  = τ0 / K
        ρ  = (n+1)/n

        # Centerline velocity
        vc  = (1/(α*ρ)) * ((abs(α)*wy - β)^ρ - (abs(α)*h - β)^ρ) - 0.5 * (dpdz / η0) * wy^2
        ∂vc = -(1/(α^2*ρ)) * ∂α * ((abs(α)*wy - β)^ρ - (abs(α)*h - β)^ρ) + 
              (1/α) * ((sign(α)*∂α*wy+abs(α)*∂wy)*(abs(α)*wy - β)^(ρ-1) - (sign(α)*∂α*h)*(abs(α)*h - β)^(ρ-1)) +
              - 0.5 * (wy^2) / η0 +
              - (dpdz / η0) * wy * ∂wy

        Q1 = (dpdz*wy^3)/(3.0*η0) + 2.0*wy*vc
        ∂Q1 = (wy^3)/(3.0*η0) + (3*dpdz*wy^2*∂wy)/(3.0*η0) + 2.0*∂wy*vc + 2.0*wy*∂vc

        Q2  = (2.0/(α*abs(α)*ρ*(ρ+1))) * ((abs(α)*h - β)^(ρ+1) - (abs(α)*wy - β)^(ρ+1)) - (2.0*(h-wy)/(ρ*α))*(abs(α)*h - β)^ρ
        ∂Q2 = -((4.0*∂α)/((α^2)*abs(α)*ρ*(ρ+1))) * ((abs(α)*h - β)^(ρ+1) - (abs(α)*wy - β)^(ρ+1)) +
              (2.0/(α*abs(α)*ρ)) * ((sign(α)*∂α*h)*(abs(α)*h - β)^(ρ) - (sign(α)*∂α*wy+abs(α)*∂wy)*(abs(α)*wy - β)^(ρ)) +
              (2.0*(∂wy)/(ρ*α))*(abs(α)*h - β)^ρ +
              (2.0*(h-wy)/(ρ*α^2))*∂α*(abs(α)*h - β)^ρ +
              - (2.0*(h-wy)/(α))*(sign(α)*∂α*h)*(abs(α)*h - β)^(ρ-1)
    end

    return Q1+Q2, ∂Q1+∂Q2 
end