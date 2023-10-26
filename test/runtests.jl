module TestingTypeUtils

using TypeUtils
using Test
using Base: OneTo
import TwoDimensional

@testset "TypeUtils" begin
    @testset "parameterless()" begin
        @test parameterless(Array) === Array
        @test parameterless(AbstractVector) === AbstractArray
        @test parameterless(DenseMatrix) === DenseArray
        @test parameterless(Array{Float32}) === Array
        @test parameterless(DenseArray{Float32,3}) === DenseArray
        @test parameterless(Matrix{Float64}) === Array
    end

    @testset "as()" begin
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
    end

    @testset "promote_eltype()" begin
        A = rand(Float32, 2, 3)
        B = ones(Int16, 2)
        C = zeros(Bool, 3, 4)
        @test promote_eltype() === promote_type()
        @test promote_eltype(A) === eltype(A)
        @test promote_eltype(B) === eltype(B)
        @test promote_eltype(C) === eltype(C)
        @test promote_eltype(A,B) === promote_type(eltype(A), eltype(B))
        @test promote_eltype(A,B,C) === promote_type(eltype(A), eltype(B), eltype(C))
    end

    @testset "convert_eltype()" begin
        # Abstract arrays.
        let A = rand(Float64, 2, 3), B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa AbstractArray{Float32,ndims(A)}
            @test B == Float32.(A)
        end
        # Base.OneTo
        let A = OneTo{Int32}(7), B = @inferred convert_eltype(Int16, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa OneTo{Int16}
            @test B == Int16.(A)
        end
        let A = OneTo(7), B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
        end
        # UnitRange
        let A = 2:8, B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
        end
        # OrdinalRange
        let A = 2.0:3.0:11.0, B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
        end
        # LinRange
        let A = LinRange(-2.0, 3.0, 5), B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa LinRange{Float32}
            @test B == Float32.(A)
        end
        # Tuples.
        let A = (1, 2, 3) #= NTuple =#, B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa NTuple{<:Any,Float32}
            @test B == Float32.(A)
        end
        let A = (1, 2.0, 3) #= not NTuple but Tuple =#, B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa NTuple{<:Any,Float32}
            @test B == Float32.(A)
        end
    end

    @testset "as_eltype()" begin
        let A = rand(Float64, 3, 4, 5), B = @inferred as_eltype(Float32, A)
            @test B == Float32.(A)
            @test eltype(B) === Float32
            @test length(B) === length(A)
            @test size(B) === size(A)
            @test axes(B) === axes(A)
            A[1,2,3] = -7
            @test B[1,2,3] === Float32(-7)
            B[1,2,3] = 19
            @test A[1,2,3] == 19
        end
        let A = view(rand(Float64, 3, 4, 5), :, 2, :), B = @inferred as_eltype(Float32, A)
            @test as_eltype(eltype(A), A) === A
            @test B == Float32.(A)
            @test eltype(B) === Float32
            @test length(B) === length(A)
            @test size(B) === size(A)
            @test axes(B) === axes(A)
            A[2,3] = -7
            @test B[2,3] === Float32(-7)
            B[2,3] = 19
            @test A[2,3] == 19
        end
        let A = 1:5
            @test as_eltype(eltype(A), A) === A
            @test as_eltype(Float32, A) == Float32.(A)
        end
        let A = 2.0:3.0:11.0
            @test as_eltype(eltype(A), A) === A
            @test as_eltype(Float32, A) == Float32.(A)
        end
    end

    @testset "return_type()" begin
        @test return_type(cos, Float32) === typeof(cos(zero(Float32)))
        @test return_type(atan, Int32) === typeof(atan(zero(Int32)))
        @test return_type(atan, Int32, Int16) === typeof(atan(zero(Int32), one(Int16)))
    end

    @testset "as_return()" begin
        f = @inferred as_return(Float32, atan)
        g = @inferred as_return(BigFloat, f)
        h = @inferred as_return(return_type(f), f)
        @test g isa parameterless(typeof(f))
        @test h === f
        @test parent(f) === atan
        @test parent(g) === parent(f)
        @test return_type(f) === Float32
        @test return_type(g) === BigFloat
        @test return_type(typeof(f)) === Float32
        @test return_type(typeof(g)) === BigFloat
        @test Base.return_types(f) === (return_type(f),)
        @test Base.return_types(g) === (return_type(g),)
        @test Base.promote_op(f, Int8) === return_type(f)
        @test Base.promote_op(f, Int8, Int16) === return_type(f)
        @test Base.promote_op(g, Float32) === return_type(g)
        @test f(0) isa Float32
        @test f(-0.3) isa Float32
        @test g(0, 1) isa BigFloat
        @test g(-0.3) isa BigFloat
    end

    # Check with TwoDimensional.
    @testset "with TwoDimensional" begin
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
end

end # module
