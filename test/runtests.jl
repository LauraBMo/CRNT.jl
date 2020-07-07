using CRNT
using Test
using Nemo

@testset "CRNT.jl" begin
    @test NonnegativeNullspaceCone([0 0; 0 0]) == [1 0; 0 1]
    @test NonnegativeNullspaceCone([1 -1; 2 -2]) == transpose([1 1])
end

@testset "TypesCompatibilites.jl" begin
    a = 0.999
    @test Float16(Nemo.fmpq(a)) == Float16(a)
    @test Nemo.fmpq(a) == 999//1000
    a = 0.9999999
    @test Float32(Nemo.fmpq(a)) == Float32(a)
    a = 0.99999999
    @test Float64(Nemo.fmpq(a)) == Float64(a)
    a = 0.11111111111111111111111
    @test_broken BigFloat(Nemo.fmpq(a)) == BigFloat(a)
    @test push!([0.1], Nemo.fmpq(1//2)) == [0.1; 0.5]
end
