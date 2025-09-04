using Test
using Documenter
using RheoJL
using LinearAlgebra
using Trapz
import Logging

@testset "Unittests" begin
    include("test_utils.jl")
    include("test_squeeze_flow.jl")
    include("test_solver.jl")
end

@testset "Doctests" begin
    # Suppress Documenter warnings during doctest
    Logging.with_logger(Logging.SimpleLogger(stderr, Logging.Error)) do
        doctest(RheoJL)
    end
end