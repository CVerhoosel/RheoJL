@testset "Sqeeze flow" begin
    @testset "Velocity profile and flux" begin
        τ₁ = 130
        K = 30
        n = 0.4
        η₀ = 1000.0
        R₀ = 0.013
        V = 0.73e-6
        h  = (V/(π*R₀^2))/2
        ∂p∂r = 4*τ₁/h
        z, v, nw = velocity_profile(h, ∂p∂r, τ₁, K, n, η₀)
        Q, ∂Q = flux(h, ∂p∂r, τ₁, K, n, η₀)

        # Trapezoidal integration flux test
        @test Q≈2*trapz(z, v) rtol=1e-6

        # Finite difference check for derivative
        ϵ     = 1e-4
        ∂p∂r₁ = (1+ϵ)*4*τ₁/h
        Q₁, _ = flux(h, ∂p∂r₁, τ₁, K, n, η₀)
        @test ∂Q*(ϵ*4*τ₁/h) ≈ (Q₁-Q) rtol=1e-3
    end

    @testset "Newtonian solution" begin
        τ₁ = Inf
        K  = 1.0
        n  = 1.0
        η₀ = 1000.0
        R₀ = 0.013
        V  = 0.73e-6
        F  = 0.5
        h₀  = (V/(π*R₀^2))/2
        N  = 100
        T  = 10.0
        γ  = 0.06
        α  = 0.5
        β  = 1e6

        # Single time step
        ∂h∂t, info = solve_system(h₀, R₀, τ₁, K, n, η₀, F, N)
        @test ∂h∂t ≈ -(8*F*h₀^3)/(3*π*η₀*R₀^4) rtol=1e-3
        @test info == 1

        # Single time step with Laplace pressure
        ∂h∂t, info = solve_system(h₀, R₀, τ₁, K, n, η₀, F, N; γ=γ, α=α)
        Fₙₜ = F - (2*γ*α/h₀)*π*R₀^2
        @test ∂h∂t ≈ -(8*Fₙₜ*h₀^3)/(3*π*η₀*R₀^4) rtol=1e-3
        @test info == 1

        # Single time step with slip
        ∂h∂t, info = solve_system(h₀, R₀, τ₁, K, n, η₀, F, N; β=β)
        ηₙₜ = η₀ / (1 + 3*η₀/(β*h₀))
        @test ∂h∂t ≈ -(8*F*h₀^3)/(3*π*ηₙₜ*R₀^4) rtol=1e-3
        @test info == 1

        # Time integration
        sol, info = integrate_system(h₀, R₀, τ₁, K, n, η₀, F, N, T; Δtₘₐₓ=0.1)
        t_end, h_end, R_end = sol[end, :]
        
        @test t_end ≈ T rtol=1e-6
        @test R_end ≈ R₀*(1 + (8*F*T*V^2)/(3*π^3*η₀*R₀^8))^(1/8) rtol=1e-4
        @test h_end ≈ (V/(π*R_end^2))/2 rtol=1e-6
        @test h_end ≈ h₀*(1 + (32*F*T*h₀^2)/(3*π*η₀*R₀^4))^(-1/4) rtol=1e-4


    end

    @testset "Force calculation" begin
        N = 100
        R₀ = 3.0
        F = 5.0

        rᵥ = collect(range(0, R₀, length=N))
        rₘ = vertex_to_midpoint(rᵥ)
        ∂p∂rₘ = -(4*F)/(π*R₀^4) * rₘ

        @test F ≈ force(rᵥ, ∂p∂rₘ) rtol=1e-3
    end
end