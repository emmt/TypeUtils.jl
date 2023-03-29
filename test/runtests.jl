module TestingAsType

using AsType
using Test
import TwoDimensional

@testset "AsType.jl" begin
    @test as(Int, 3) === 3
    @test as(Int, 1.0) === 1
    @test as(Int16, 1.0) === Int16(1)
    @test_throws InexactError as(Int, sqrt(2))
    @test_throws Exception as(Int, π)
    @test map(as(Int), (Int8(1), Int16(2), Int32(3), Int64(4))) === (1,2,3,4)

    @test as(Tuple, CartesianIndex((1,2,3))) === (1,2,3)
    @test as(NTuple{3}, CartesianIndex((1,2,3))) === (1,2,3)
    @test as(Tuple, CartesianIndices((2:5, -1:4))) === (2:5, -1:4)
    @test as(NTuple{2}, CartesianIndices((2:5, -1:4))) === (2:5, -1:4)

    @test as(CartesianIndex,()) === CartesianIndex()
    @test as(CartesianIndex, CartesianIndex((1,2,3))) === CartesianIndex((1,2,3))
    @test as(CartesianIndex{3}, CartesianIndex((1,2,3))) === CartesianIndex((1,2,3))
    @test as(CartesianIndex, (0x1,2,Int16(3))) === CartesianIndex((1,2,3))
    @test as(CartesianIndex{3}, (0x1,2,Int16(3))) === CartesianIndex((1,2,3))

    @test as(CartesianIndices,()) === CartesianIndices(())
    @test as(CartesianIndices,(2:3,6)) === CartesianIndices((2:3,6))
    @test as(CartesianIndices,CartesianIndices((2:3,6))) === CartesianIndices((2:3,6))
    @test as(CartesianIndices{3},(2:3,6,-1:4)) === CartesianIndices((2:3,6,-1:4))
    @test as(CartesianIndices{3},CartesianIndices((2:3,6,-1:4))) === CartesianIndices((2:3,6,-1:4))

    @test as(String, :hello) == "hello"
    @test as(String, :hello) isa String
    @test as(Symbol, "hello") === :hello

    # Check with TwoDimensional.
    let Point = TwoDimensional.Point,
        WeightedPoint = TwoDimensional.WeightedPoint,
        BoundingBox = TwoDimensional.BoundingBox

        @test as(Tuple, Point(11, -9)) === (11, -9)
        @test as(NTuple{2,Int16}, Point(11, -9)) === (Int16(11), Int16(-9))
        @test as(Point, (11, -9)) === Point(11, -9)
        @test as(Point{Float32}, (11, -9)) === Point{Float32}(11, -9)

        @test as(Tuple, WeightedPoint(2.0, 11.0, -9.0)) === (2.0, 11.0, -9.0)
        @test as(NTuple{3,Float32}, WeightedPoint(2.0, 11.0, -9.0)) === (Float32(2), Float32(11), Float32(-9))
        @test as(WeightedPoint, (2.0, 11.0, -9.0)) === WeightedPoint(2.0, 11.0, -9.0)
        @test as(WeightedPoint{Float32}, (2.0, 11.0, -9.0)) === WeightedPoint{Float32}(2.0, 11.0, -9.0)

        @test as(Tuple, BoundingBox(2, 11, -9, 7)) === (2, 11, -9, 7)
        @test as(NTuple{4,Int16}, BoundingBox(2, 11, -9, 7)) === map(Int16, (2, 11, -9, 7))
        @test as(BoundingBox, (2, 11, -9, 7)) === BoundingBox(2, 11, -9, 7)
        @test as(BoundingBox{Float32}, (2, 11, -9, 7)) === BoundingBox{Float32}(2, 11, -9, 7)
    end
end

end # module
