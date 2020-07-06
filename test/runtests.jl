using CRNT
using Test

@testset "CRNT.jl" begin
    @test NonnegativeNullspaceCone([0 0; 0 0]) == [1 0; 0 1]
    @test NonnegativeNullspaceCone([1 -1; 2 -2]) == transpose([1 1])
end
