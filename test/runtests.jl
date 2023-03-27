module TestingAsType

using AsType
using Test

@testset "AsType.jl" begin
    @test as(Int, 3) === 3
    @test as(Int, 1.0) === 1
    @test as(Int16, 1.0) === Int16(1)
    @test_throws InexactError as(Int, sqrt(2))
    @test_throws Exception as(Int, π)
    @test map(as(Int), (Int8(1), Int16(2), Int32(3), Int64(4))) === (1,2,3,4)

    @test as(CartesianIndex,()) === CartesianIndex()
    @test as(CartesianIndex, CartesianIndex((1,2,3))) === CartesianIndex((1,2,3))
    @test as(CartesianIndex{3}, CartesianIndex((1,2,3))) === CartesianIndex((1,2,3))
    @test as(CartesianIndex, (0x1,2,Int16(3))) === CartesianIndex((1,2,3))
    @test as(CartesianIndex{3}, (0x1,2,Int16(3))) === CartesianIndex((1,2,3))

    @test as(Tuple, CartesianIndex((1,2,3))) === (1,2,3)
    @test as(Tuple, CartesianIndices((2:5, -1:4))) === (2:5, -1:4)
end

end # module
