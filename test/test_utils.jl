@testset "Integration utilities" begin
    N = 100
    R₀ = 3.0
    F = 5.0

    rᵥ = collect(range(0, R₀, length=N))
    rₘ = vertex_to_midpoint(rᵥ)
    pᵥ₀ = (2*F)/(π*R₀^2) * (1.0 .- (rᵥ ./ R₀).^2)
    ∂p∂rₘ₀ = -(4*F)/(π*R₀^4) * rₘ

    @test norm(pᵥ₀ .- right_integrate(rᵥ, ∂p∂rₘ₀)) ≈ 0.0 atol=1e-10
    @test integral(rᵥ, vertex_to_midpoint(2*π*rᵥ .* pᵥ₀)) ≈ F rtol=1e-3 
end