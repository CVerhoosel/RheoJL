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
    @test Q≈2*trapz(z, v) rtol=1e-6
end