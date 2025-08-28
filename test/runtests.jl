using Test
using Documenter
using RheoJL

@testset "Unittests" begin
    @test vertex_to_midpoint([1.0, 3.0, 7.0]) == [2.0, 5.0]
end

@testset "Doctests" begin
    doctest(RheoJL)
end