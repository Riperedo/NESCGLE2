# test/runtests.jl
using NESCGLE2
using Test

@testset "ThermodynamicInstabilities" begin
    # Tus pruebas para las funciones de este módulo
    # Ejemplo:
    # result = NESCGLE2.calcular_spinodal(...)
    # @test result ≈ expected_value atol=1e-6
    include("testThermodynamicInstabilities.jl")
end

@testset "AsymptoticProperties" begin
    # ...
end

@testset "StructuralRelaxation" begin
    # ...
end