@testset "Sqeeze flow" begin
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