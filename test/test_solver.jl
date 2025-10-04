@testset "Newton solver" begin
    
    @testset "Linear function" begin
        r = x -> (3x - 5, 3, [])

        x, info = newton(r, x₀=0.0)
        @test x ≈ 5/3
        @test info == 1

        x, info = newton(r, x₀=5/3)
        @test x ≈ 5/3
        @test info == 0

        x, info = newton(r, x₀=0.0, bracket=[0.0, 5.0])
        @test x ≈ 5/3
        @test info == 1

        x, info = newton(r, x₀=5, bracket=[5/3, 5.0])
        @test x ≈ 5/3
        @test info == 0

        x, info = newton(r, x₀=0.0, bracket=[0.0, 5/3])
        @test x ≈ 5/3
        @test info == 0
    end
   
   @testset "Nonlinear function" begin 
        r = x -> (x^3 - 2x - 5, 3x^2 - 2, [])

        x, info = newton(r, x₀=4.0, bracket=[0.0,5.0])
        @test x ≈ 2.094551525517854 rtol=1e-7
        @test info == 6
   end
end