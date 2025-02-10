module TestingTypeUtilsWithoutExtensions
using TypeUtils
using Test

if isempty(filter(x -> x.name == "OffsetArrays", keys(Base.loaded_modules)))
    @testset "TypeUtils (without extensions)" begin
        @testset "new_array()" begin
            @test_throws Exception new_array(Int16, -1:2)
        end
    end
end

end # module TestingTypeUtilsWithoutExtensions

module TestingTypeUtils

using TypeUtils
using TypeUtils: BareNumber
using Unitful
using OffsetArrays
using Test
using Base: OneTo

struct Foo{T1,T2}
    z::Complex{T1}
    r::T2
    i::Int
end

struct Bar{T}
    x::T
    y::Tuple{T,Int16,UInt8}
end

# Minimal implementation of a custom numeric type.
struct MyNumber{T<:BareNumber} <: Number
    val::T
end
Base.zero(::Type{MyNumber{T}}) where {T} = MyNumber{T}(zero(T))
Base.oneunit(::Type{MyNumber{T}}) where {T} = MyNumber{T}(one(T))
Base.one(::Type{MyNumber{T}}) where {T} = one(T)
Base.real(x::MyNumber{<:Real}) = x
Base.real(x::MyNumber{<:Complex}) = MyNumber(real(x.val))
Base.real(::Type{MyNumber{T}}) where {T<:Real} = MyNumber{T}
Base.real(::Type{MyNumber{Complex{T}}}) where {T<:Real} = MyNumber{T}
TypeUtils.unitless(x::MyNumber) = x.val
TypeUtils.convert_bare_type(T::Type{<:BareNumber}, ::Type{<:MyNumber}) =
    MyNumber{T}
TypeUtils.convert_real_type(T::Type{<:Real}, ::Type{MyNumber{S}}) where {S} =
    MyNumber{convert_real_type(T,S)}

# Different implementations of in-place multiplication.
function scale!(::Val{1}, A::AbstractArray, α::Number)
    alpha = convert_bare_type(eltype(A), α)
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end

function scale!(::Val{2}, A::AbstractArray, α::Number)
    alpha = convert_real_type(eltype(A), α)
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end

function scale!(::Val{3}, A::AbstractArray, α::Union{Real,Complex})
    alpha = if α isa Complex && bare_type(eltype(A)) isa Real
        convert(real_type(eltype(A)), α)
    else
        convert_real_type(eltype(A), α)
    end
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end

# Type-instable function.
type_instable_function(i::Integer) = isodd(i) ? π : sqrt(i)

@testset "TypeUtils" begin
    @testset "Miscellaneous" begin
        # Check that TypeUtils.Unsupported cannot be instantiated.
        @test_throws Exception TypeUtils.Unsupported()
    end

    @testset "parameterless()" begin
        @test parameterless(Array) === Array
        @test parameterless(AbstractVector) === AbstractArray
        @test parameterless(DenseMatrix) === DenseArray
        @test parameterless(Array{Float32}) === Array
        @test parameterless(DenseArray{Float32,3}) === DenseArray
        @test parameterless(Matrix{Float64}) === Array
        @test parameterless(typeof((1, 2, 3))) === Tuple
        @test parameterless(typeof((1, π, "e"))) === Tuple
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

    @testset "Array axes and size" begin
        @test TypeUtils.Dim === Int
        @test  (42                   isa eltype(ArrayShape))
        @test  (42                   isa eltype(RelaxedArrayShape))
        @test  (Int8(42)             isa eltype(ArrayShape))
        @test  (Int8(42)             isa eltype(RelaxedArrayShape))
        @test  (-1:5                 isa eltype(ArrayShape))
        @test  (-1:5                 isa eltype(RelaxedArrayShape))
        @test  (0x2:0x5              isa eltype(ArrayShape))
        @test  (0x2:0x5              isa eltype(RelaxedArrayShape))
        @test !(-1:1:5               isa eltype(ArrayShape))
        @test  (-1:1:5               isa eltype(RelaxedArrayShape))
        @test  (Base.OneTo(9)        isa eltype(ArrayShape))
        @test  (Base.OneTo(9)        isa eltype(RelaxedArrayShape))
        @test  (Base.OneTo{Int16}(9) isa eltype(ArrayShape))
        @test  (Base.OneTo{Int16}(9) isa eltype(RelaxedArrayShape))

        @test_throws MethodError as_array_dim(π)
        @test as_array_dim(7) === 7
        @test as_array_dim(Int16(9)) === 9
        @test as_array_dim(Base.OneTo(5)) === 5
        @test as_array_dim(-1:4) === 6
        @test as_array_dim(Int8(0):Int8(3)) === 4
        @test as_array_dim(-1:1:4) === 6
        @test_throws ArgumentError as_array_dim(-1:2:3) # non-unit step

        @test_throws MethodError as_array_axis(π)
        @test as_array_axis(7) === Base.OneTo(7)
        @test as_array_axis(Int16(9)) === Base.OneTo(9)
        @test as_array_axis(Base.OneTo(5)) === Base.OneTo(5)
        @test as_array_axis(-1:4) === -1:4
        @test as_array_axis(1:4) === 1:4
        @test as_array_axis(Int8(0):Int8(3)) === 0:3
        @test as_array_axis(-1:1:4) === -1:4
        @test_throws ArgumentError as_array_axis(-1:2:3) # non-unit step

        @test_throws MethodError as_array_shape(π)
        @test_throws MethodError as_array_shape((π,))
        @test as_array_shape() === ()
        @test as_array_shape(()) === ()
        @test as_array_shape(5) === (5,)
        @test as_array_shape(2:5) === (2:5,)
        @test as_array_shape(2,Int16(9),Int8(3)) === (2,9,3)
        @test as_array_shape((2,Int16(9),Int8(3))) === (2,9,3)
        @test as_array_shape(2,3,4) === (2,3,4)
        @test as_array_shape(2,Base.OneTo{Int8}(3),4) === (2,3,4)
        @test as_array_shape(2,Int8(1):Int8(3),4) === (Base.OneTo(2),1:3,Base.OneTo(4))
        @test as_array_shape((0:2,-4:4,-2:1)) === (0:2,-4:4,-2:1)
        @test_throws ArgumentError as_array_shape((1:4, 1:2:6,)) # non-unit step

        @test_throws MethodError as_array_axes(π)
        @test_throws MethodError as_array_axes((π,))
        @test as_array_axes() === ()
        @test as_array_axes(()) === ()
        @test as_array_axes(5) === (Base.OneTo(5),)
        @test as_array_axes(2:5) === (2:5,)
        @test as_array_axes(2,Int16(9),Int8(3)) === (Base.OneTo(2),Base.OneTo(9),Base.OneTo(3))
        @test as_array_axes((2,Int16(9),Int8(3))) === (Base.OneTo(2),Base.OneTo(9),Base.OneTo(3))
        @test as_array_axes(2,3,4) === (Base.OneTo(2),Base.OneTo(3),Base.OneTo(4))
        @test as_array_axes(2,Base.OneTo{Int8}(3),4) === (Base.OneTo(2),Base.OneTo(3),Base.OneTo(4))
        @test as_array_axes(2,Int8(1):Int8(3),4) === (Base.OneTo(2),1:3,Base.OneTo(4))
        @test as_array_axes((0:2,-4:4,-2:1)) === (0:2,-4:4,-2:1)
        @test_throws ArgumentError as_array_axes((1:4, 1:2:6,)) # non-unit step

        @test_throws MethodError as_array_size(π)
        @test_throws MethodError as_array_size((π,))
        @test as_array_size() === ()
        @test as_array_size(()) === ()
        @test as_array_size(5) === (5,)
        @test as_array_size(2:5) === (4,)
        @test as_array_size(2,Int16(9),Int8(3)) === (2,9,3)
        @test as_array_size((2,Int16(9),Int8(3))) === (2,9,3)
        @test as_array_size(2,3,4) === (2,3,4)
        @test as_array_size(2,Base.OneTo{Int8}(3),4) === (2,3,4)
        @test as_array_size(2,Int8(1):Int8(3),4) === (2,3,4)
        @test as_array_size((0:2,-4:4,-2:1)) === (3,9,4)
        @test_throws ArgumentError as_array_size((1:4, 1:2:6,)) # non-unit step
    end

    @testset "new_array()" begin
        A = @inferred new_array(Float32, 2,3,4)
        @test A isa Array{Float32,3}
        @test size(A) === (2,3,4)
        B = @inferred new_array(Float64, axes(A))
        @test B isa Array{Float64,3}
        @test size(B) === size(A)
        C = @inferred new_array(Int16, 2,-1:1,4)
        @test C isa OffsetArray{Int16,3}
        @test size(C) === (2,3,4)
        @test axes(C) == (1:2, -1:1, 1:4)
        D = @inferred new_array(Int8, axes(C))
        @test D isa OffsetArray{Int8,3}
        @test axes(D) == axes(C)
    end

    @testset "to_same_concrete_type()" begin
        @test to_same_concrete_type(Int) === Int
        @test to_same_concrete_type(UInt8, UInt8) === UInt8
        @test to_same_concrete_type(Int8, Int8, Int8) === Int8
        @test to_same_concrete_type(UInt8, UInt16) === UInt16
        @test to_same_concrete_type(Int8, Int16, Int32) === Int32
        @test_throws ArgumentError to_same_concrete_type()
        @test_throws ArgumentError to_same_concrete_type(AbstractUnitRange)
        @test_throws ArgumentError to_same_concrete_type(AbstractFloat, Float32)
        @test_throws ArgumentError to_same_concrete_type(String, Float32, Int16)
    end

    @testset "to_same_type()" begin
        @test to_same_type() === ()
        @test to_same_type(pi) === (pi,)
        @test to_same_type('x') === ('x',)
        @test to_same_type(1, 2) === (1, 2)
        @test to_same_type(0x1, 0x2, 0x3) === (0x1, 0x2, 0x3)
        @test to_same_type(0x1, 2) === (1, 2)
        @test to_same_type(UInt8(1), Int8(2), Int16(3)) === (Int16(1), Int16(2), Int16(3))
        @test_throws ArgumentError to_same_type(2u"mm", 3.0)
        @test_throws ArgumentError to_same_type(2u"mm", 3.0u"kg")
        x1, x2 = @inferred to_same_type(2u"mm", 3.0u"cm")
        @test typeof(x1) === typeof(x2)
        @test_throws ArgumentError to_same_type(Int)
        @test_throws ArgumentError to_same_type(Int8,Int8)
        @test_throws ArgumentError to_same_type(Int8,Int16)
        @test_throws ArgumentError to_same_type(Int16,Int16,Int16)
        @test_throws ArgumentError to_same_type(Int16,Int32,Int64)
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
        # Numbers.
        let A = rand(Float64), B = @inferred convert_eltype(Float32, A)
            @test convert_eltype(eltype(A), A) === A
            @test B isa Float32
            @test B == Float32.(A)
        end
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
        @test convert(TypeStableFunction, f) === f
        @test convert(TypeStableFunction{return_type(f)}, f) === f
        @test convert(TypeStableFunction{return_type(g)}, f) === g
        @test convert(AbstractTypeStableFunction, f) === f
        @test convert(AbstractTypeStableFunction{return_type(f)}, f) === f
        @test convert(AbstractTypeStableFunction{return_type(g)}, f) === g

        # Make a type-stable function out of a type-instable one.
        @test typeof(type_instable_function(1)) != typeof(type_instable_function(2))
        if VERSION < v"1.1"
            # `Base.promote_op` returns `Any` for Julia version prior to 1.1 so we must
            # provide the return type ourself.
            f = @inferred TypeStableFunction{Float64}(type_instable_function)
        else
            f = @inferred TypeStableFunction(type_instable_function, Int)
            @test AbstractTypeStableFunction(type_instable_function, Int) === f
        end
        @test parent(f) === type_instable_function
        @test return_type(typeof(f)) === Float64
        @test typeof(f(1)) === typeof(f(2)) === return_type(typeof(f))
        @test TypeStableFunction{return_type(typeof(f))}(type_instable_function) === f
        @test AbstractTypeStableFunction{return_type(typeof(f))}(type_instable_function) === f
    end

    @testset "Numeric types" begin
        # bare_type with no argument
        @test bare_type() === BareNumber

        # bare_type for values
        @test bare_type(1.0) === Float64
        @test bare_type(Float32) === Float32
        @test bare_type(Complex(2,3)) === Complex{Int}
        @test bare_type(NaN) === typeof(NaN)
        @test bare_type(π) === typeof(π)
        @test bare_type(3//4) === typeof(3//4)
        @test_throws ErrorException bare_type("hello")

        # bare_type for types
        @test bare_type(Real) === Real
        @test bare_type(Integer) === Integer
        @test bare_type(Float32) === Float32
        @test bare_type(BigFloat) === BigFloat
        @test bare_type(Complex{Int}) === Complex{Int}
        @test bare_type(typeof(π)) === typeof(π)
        @test bare_type(typeof(3//4)) === typeof(3//4)
        @test_throws ErrorException bare_type(AbstractString)

        # bare_type with multiple arguments
        @test bare_type(1, 0f0) === Float32
        @test bare_type(Int, pi) === promote_type(Int, typeof(pi))
        @test bare_type(4, pi, 1.0) === promote_type(Int, typeof(pi), Float64)
        @test bare_type(Int, Int8, Float32) === promote_type(Int, Int8, Float32)
        @test bare_type(Int, Int8, Float32) === promote_type(Int, Int8, Float32)
        @test bare_type(Int, Int8, Int16, Float32) === promote_type(Int, Int8, Int16, Float32)

        # default implementation
        @test bare_type(MyNumber(1.2f0)) === Float32
        @test bare_type(MyNumber{Int16}) === Int16

        # convert_bare_type
        @test convert_bare_type(Int, -1) === -1
        @test convert_bare_type(Int, 2.0) === 2
        @test convert_bare_type(Float32, 2.0) === 2.0f0
        @test convert_bare_type(MyNumber{Int16}, 12.0) === Int16(12)
        @test_throws ErrorException convert_bare_type(Int, "oups!")
        for x in (missing, nothing, undef)
            @test convert_bare_type(Int8, x) === x
            @test convert_bare_type(Int8, typeof(x)) === typeof(x)
        end
        @test convert_bare_type(Int16, Int32) === Int16
        @test convert_bare_type(Int16, Complex{Int32}) === Int16
        @test convert_bare_type(Complex{Int16}, Int32) === Complex{Int16}
        @test convert_bare_type(Complex{Int16}, Complex{Int32}) === Complex{Int16}
        @test convert_bare_type(MyNumber{Int16}, MyNumber{Int32}) === MyNumber{Int16}
        @test convert_bare_type(Int16, MyNumber{Complex{Int32}}) === MyNumber{Int16}
        @test convert_bare_type(Complex{Int16}, MyNumber{Int32}) === MyNumber{Complex{Int16}}
        @test_throws ErrorException convert_bare_type(Int, String)

        # real_type with no argument
        @test real_type() === Real

        # real_type for values
        @test real_type(1.0) === Float64
        @test real_type(Float32) === Float32
        @test real_type(Complex(2,3)) === Int
        @test real_type(NaN) === typeof(NaN)
        @test real_type(π) === typeof(π)
        @test real_type(3//4) === typeof(3//4)
        @test_throws ErrorException real_type("hello")

        # real_type for types
        @test real_type(Real) === Real
        @test real_type(Integer) === Integer
        @test real_type(Float32) === Float32
        @test real_type(BigFloat) === BigFloat
        @test real_type(Complex{Int}) === Int
        @test real_type(typeof(π)) === typeof(π)
        @test real_type(typeof(3//4)) === typeof(3//4)
        @test_throws ErrorException real_type(AbstractString)

        # real_type with multiple arguments
        @test real_type(1, 0f0) === Float32
        @test real_type(Int, pi) === promote_type(Int, typeof(pi))
        @test real_type(4, pi, 1.0) === promote_type(Int, typeof(pi), Float64)
        @test real_type(Int, Int8, Float32) === promote_type(Int, Int8, Float32)
        @test real_type(Int, Int8, Float32) === promote_type(Int, Int8, Float32)
        @test real_type(Int, Int8, Complex{Int16}, Float32) === promote_type(Int, Int8, Int16, Float32)

        # default implementation
        @test real_type(MyNumber(1.2f0)) === Float32
        @test real_type(MyNumber{Complex{Int16}}) === Int16

        # convert_real_type
        @test convert_real_type(Int, -1) === -1
        @test convert_real_type(Int, 2.0) === 2
        @test convert_real_type(Complex{Float32}, 2.0) === 2.0f0
        let x = 2.1f0
            @test convert_real_type(typeof(x), x) === x
        end
        let x = 2 + 3im
            @test convert_real_type(typeof(x), x) === x
            @test convert_real_type(real(typeof(x)), x) === x
        end
        @test convert_real_type(Complex{Int16}, 2 + 3im) === Complex{Int16}(2, 3)
        @test convert_real_type(Float32, 2.0 - 1.0im) === Complex{Float32}(2, -1)
        @test convert_real_type(MyNumber{Int16}, 12.0) === Int16(12)
        @test_throws ErrorException convert_real_type(Int, "oups!")
        for x in (missing, nothing, undef)
            @test convert_real_type(Int8, x) === x
            @test convert_real_type(Int8, typeof(x)) === typeof(x)
        end
        @test convert_real_type(Int16, Int32) === Int16
        @test convert_real_type(Int16, Complex{Int32}) === Complex{Int16}
        @test convert_real_type(Complex{Int16}, Int32) === Int16
        @test convert_real_type(Complex{Int16}, Complex{Int32}) === Complex{Int16}
        @test convert_real_type(MyNumber{Int16}, MyNumber{Int32}) === MyNumber{Int16}
        @test convert_real_type(Int16, MyNumber{Complex{Int32}}) === MyNumber{Complex{Int16}}
        @test convert_real_type(Complex{Int16}, MyNumber{Int32}) === MyNumber{Int16}
        @test_throws ErrorException convert_real_type(Int, String)

        # floating_point_type
        @test floating_point_type() === AbstractFloat
        @test floating_point_type(Int16) === float(Int16)
        @test floating_point_type(Int16, Float32) === Float32
        @test floating_point_type(Int16, Float32, Complex{Float64}) === Float64
        @test floating_point_type(Int16, Float32, Complex{Float16}, 0x1) === Float32

        # convert_floating_point_type
        @test convert_floating_point_type(Int, -1) === -1.0
        @test convert_floating_point_type(Int, 2.0) === 2.0
        @test convert_floating_point_type(Complex{Float32}, 2.0) === 2.0f0
        let x = 2.1f0
            @test convert_floating_point_type(typeof(x), x) === x
        end
        let x = 2 + 3im
            @test convert_floating_point_type(typeof(x), x) === float(x)
            @test convert_floating_point_type(real(typeof(x)), x) === float(x)
        end
        @test convert_floating_point_type(Complex{Int16}, 2 + 3im) === Complex{Float64}(2, 3)
        @test convert_floating_point_type(Float32, 2.0 - 1.0im) === Complex{Float32}(2, -1)
        @test convert_floating_point_type(MyNumber{Int16}, 12.0) === Float64(12)
        @test_throws ErrorException convert_floating_point_type(Int, "oups!")
        for x in (missing, nothing, undef)
            @test convert_floating_point_type(Float32, x) === x
            @test convert_floating_point_type(Float32, typeof(x)) === typeof(x)
        end
        @test convert_floating_point_type(Float32, Float64) === Float32
        @test convert_floating_point_type(Float32, Complex{Float64}) === Complex{Float32}
        @test convert_floating_point_type(Complex{Float32}, Float64) === Float32
        @test convert_floating_point_type(Complex{Float32}, Complex{Float64}) === Complex{Float32}
        @test convert_floating_point_type(MyNumber{Float32}, MyNumber{Float64}) === MyNumber{Float32}
        @test convert_floating_point_type(Float32, MyNumber{Complex{Float64}}) === MyNumber{Complex{Float32}}
        @test convert_floating_point_type(Complex{Float32}, MyNumber{Float64}) === MyNumber{Float32}
        @test_throws ErrorException convert_floating_point_type(Int, String)

        # unitless
        @test unitless(Real) === bare_type(Real)
        @test unitless(Integer) === bare_type(Integer)
        @test unitless(Float32) === bare_type(Float32)
        @test unitless(BigFloat) === bare_type(BigFloat)
        @test unitless(Complex{Int}) === bare_type(Complex{Int})
        @test unitless(typeof(3//4)) === bare_type(typeof(3//4))
        @test unitless(typeof(π)) === bare_type(typeof(π))
        @test unitless(17.0) === 17.0
        @test unitless(17.0f0) === 17.0f0
        @test unitless(17) === 17
        @test unitless(Int16(17)) === Int16(17)
        @test unitless(true) === true
        @test unitless(false) === false
        @test unitless(3//4) === 3//4
        @test unitless(π) === π
        @test unitless(MyNumber(1.2f0)) === 1.2f0
        @test unitless(MyNumber{Int16}) === Int16

        # in-place multiplication
        A = [1.1, 1.3, 2.7]
        B = similar(A)
        alpha = 2 + 0im
        @test scale!(Val(1), copyto!(B, A), alpha) == alpha .* A
        @test scale!(Val(2), copyto!(B, A), alpha) == alpha .* A
        @test scale!(Val(3), copyto!(B, A), alpha) == alpha .* A
        alpha = 2 - 1im
        @test_throws InexactError scale!(Val(1), copyto!(B, A), alpha)
        @test_throws InexactError scale!(Val(2), copyto!(B, A), alpha)
        @test_throws InexactError scale!(Val(3), copyto!(B, A), alpha)

    end

    @testset "Unitful quantities" begin
        # bare_type for values
        @test bare_type(u"2.0m/s") === Float64
        @test bare_type(u"35GHz") === Int

        # bare_type for types
        @test bare_type(typeof(u"2.0m/s")) === Float64
        @test bare_type(typeof(u"35GHz")) === Int

        # bare_type with multiple arguments
        @test bare_type(u"2.0m/s", u"35GHz") === Float64
        @test bare_type(1, u"2.0f0m/s", u"35GHz") === Float32
        @test bare_type(1, u"2.0f0m/s", u"35GHz", Complex{Int8}(11)) === Complex{Float32}

        # convert_bare_type
        @test convert_bare_type(Float64, u"2.0m/s") === u"2.0m/s"
        @test convert_bare_type(Int, u"2.0m/s") === u"2m/s"
        @test convert_bare_type(Float32, u"35GHz") === u"35.0f0GHz"
        let T = typeof(u"3.5GHz")
            for x in (missing, nothing, undef)
                @test convert_bare_type(T, x) === x
                @test convert_bare_type(T, typeof(x)) === typeof(x)
            end
        end
        let u = u"km/s"
            @test convert_bare_type(Int16, typeof(one(Int32)*u)) === typeof(one(Int16)*u)
            @test convert_bare_type(Int16, typeof(one(Complex{Int32})*u)) === typeof(one(Int16)*u)
            @test convert_bare_type(Complex{Int16}, typeof(one(Int32)*u)) === typeof(one(Complex{Int16})*u)
            @test convert_bare_type(Complex{Int16}, typeof(one(Complex{Int32})*u)) === typeof(one(Complex{Int16})*u)
        end

        # real_type for values
        @test real_type(u"2.0m/s") === Float64
        @test real_type(u"35GHz") === Int

        # real_type for types
        @test real_type(typeof(u"2.0m/s")) === Float64
        @test real_type(typeof(u"35GHz")) === Int

        # real_type with multiple arguments
        @test real_type(u"2.0m/s", u"35GHz") === Float64
        @test real_type(1, u"2.0f0m/s", u"35GHz") === Float32
        @test real_type(1, u"2.0f0m/s", u"35GHz", Complex{Int8}(11)) === Float32

        # convert_real_type
        @test convert_real_type(Float64, u"2.0m/s") === u"2.0m/s"
        @test convert_real_type(Int, u"2.0m/s") === u"2m/s"
        @test convert_real_type(Float32, u"35GHz") === u"35.0f0GHz"
        let T = typeof(u"3.5GHz")
            for x in (missing, nothing, undef)
                @test convert_real_type(T, x) === x
                @test convert_real_type(T, typeof(x)) === typeof(x)
            end
        end
        let u = u"km/s"
            @test convert_real_type(Int16, typeof(one(Int32)*u)) === typeof(one(Int16)*u)
            @test convert_real_type(Int16, typeof(one(Complex{Int32})*u)) === typeof(one(Complex{Int16})*u)
            @test convert_real_type(Complex{Int16}, typeof(one(Int32)*u)) === typeof(one(Int16)*u)
            @test convert_real_type(Complex{Int16}, typeof(one(Complex{Int32})*u)) === typeof(one(Complex{Int16})*u)
        end

        # convert_floating_point_type
        @test convert_floating_point_type(Float64, u"2.0m/s") === u"2.0m/s"
        @test convert_floating_point_type(Int, u"2.0m/s") === u"2.0m/s"
        @test convert_floating_point_type(Float32, u"35GHz") === u"35.0f0GHz"
        let T = typeof(u"3.5GHz")
            for x in (missing, nothing, undef)
                @test convert_floating_point_type(T, x) === x
                @test convert_floating_point_type(T, typeof(x)) === typeof(x)
            end
        end
        let u = u"km/s"
            @test convert_floating_point_type(Int16, typeof(one(Int32)*u)) === typeof(float(one(Int16))*u)
            @test convert_floating_point_type(Int16, typeof(one(Complex{Int32})*u)) === typeof(float(one(Complex{Int16}))*u)
            @test convert_floating_point_type(Complex{Int16}, typeof(one(Int32)*u)) === typeof(float(one(Int16))*u)
            @test convert_floating_point_type(Complex{Int16}, typeof(one(Complex{Int32})*u)) === typeof(float(one(Complex{Int16}))*u)
        end

        # unitless
        @test unitless(u"17GHz") === 17
        @test unitless(typeof(u"2.0f0m/s")) === Float32

        # in-place multiplication
        A = [1.1u"m/s", 1.3u"m/s", 2.7u"m/s"]
        B = similar(A)
        alpha = 2 + 0im
        @test scale!(Val(1), copyto!(B, A), alpha) == alpha .* A
        @test scale!(Val(2), copyto!(B, A), alpha) == alpha .* A
        @test scale!(Val(3), copyto!(B, A), alpha) == alpha .* A
        alpha = 2 - 1im
        @test_throws InexactError scale!(Val(1), copyto!(B, A), alpha)
        @test_throws InexactError scale!(Val(2), copyto!(B, A), alpha)
        @test_throws InexactError scale!(Val(3), copyto!(B, A), alpha)
    end

    @testset "Destructure and restructure" begin
        # Test with structures.
        obj = Foo(2.0f0 - 3.0f0im, π, 42) # NOTE π is special
        @test obj isa Foo{Float32,typeof(π)}
        @test isconcretetype(typeof(obj))
        @test struct_length(obj) == 4
        vals = @inferred destructure(obj)
        @test vals === (2.0f0, -3.0f0, π, 42)
        @test obj === @inferred restructure(typeof(obj), vals)
        @test obj === restructure(parameterless(typeof(obj)), vals)
        @test obj === @inferred restructure(typeof(obj), (0, 1, vals...); offset=2)
        vec = Vector{Float64}(undef, 6)
        @test destructure!(vec, obj)[1:4] ≈ collect(vals)
        @test destructure!(vec, obj; offset=1)[2:5] ≈ collect(vals)
        x = restructure(Foo{Float64,Float64}, vec; offset=1)
        @test x isa Foo{Float64,Float64}
        @test x.z ≈ obj.z
        @test x.r ≈ obj.r
        @test x.i == obj.i
        # Test with tuples.
        a = (1, 2, (3, (4, 5)), π, ())
        @test struct_length(a) == 7
        let vals = @inferred destructure(a)
            if VERSION < v"1.0" || v"1.5" < VERSION < v"1.10"
                # FIXME: not all versions of Julia corectly infer the result.
                @test a === restructure(typeof(a), vals)
            else
                @test a === @inferred restructure(typeof(a), vals)
            end
        end
        # Test with structure and tuples.
        b = (1, 2, obj, (3, 4))
        @test struct_length(b) == 4 + struct_length(obj)
        let vals = @inferred destructure(b)
            if VERSION < v"1.0" || v"1.5" < VERSION < v"1.10"
                # FIXME: not all versions of Julia corectly infer the result.
                @test b === restructure(typeof(b), vals)
            else
                @test b === @inferred restructure(typeof(b), vals)
            end
        end
        c = Bar{Float32}(1, (2, 3, 4))
        @test struct_length(c) == 4
        let vals = @inferred destructure(c)
            @test c === @inferred restructure(typeof(c), vals)
        end
    end
end

end # module TestingTypeUtils
