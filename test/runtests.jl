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
using TypeUtils: BareNumber, BIT_INTEGERS, Unsupported
using LinearAlgebra
using Unitful
using OffsetArrays
using Test
using Base: OneTo

"""
    x ≗ y

yields whether `x` and `y` have the same element types, the same axes, and the same values
(in the sense of `isequal`). It can be seen as a shortcut for:

    eltype(x) == eltype(y) && axes(x) == axes(y) && all(isequal, x, y)

"""
≗(x::Any, y::Any) = false
≗(x::T, y::T) where {T} = isequal(x, y)
function ≗(x::AbstractArray{T,N}, y::AbstractArray{T,N}) where {T,N}
    axes(x) == axes(y) || return false
    @inbounds for i in eachindex(x, y)
        isequal(x[i], y[i]) || return false
    end
    return true
end

# Check for approximate equality for numbers, arrays, and factorizations.
≃(x::Any, y::Any; kwds...) = isequal(x, y)
for type in (Number, AbstractArray, LinearAlgebra.AbstractQ,)
    @eval ≃(x::$type, y::$type; kwds...) = isapprox(x, y; kwds...)
end
function ≃(x::Factorization, y::Factorization; kwds...)
    x === y && return true
    parameterless(typeof(x)) == parameterless(typeof(y)) || return false
    axes(x) == axes(y) || return false
    # NOTE Must compare properties, not fields and do not use `hasproperty` which only
    # appears in Julia 1.2.
    ppt_x = propertynames(x)
    ppt_y = propertynames(y)
    for key in ppt_x
        key ∈ ppt_y || return false
    end
    for key in ppt_y
        key ∈ ppt_x || return false
        ≃(getproperty(x, key), getproperty(y, key); kwds...) || return false
    end
    return true
end

struct TestUnitRange{T<:Real} <: AbstractUnitRange{T}
    length::Int
    start::T
    function TestUnitRange{T}(; start::Real, length::Integer) where {T<:Real}
        length ≥ 0 || error("length must be ≥ 0")
        return new{T}(length, start)
    end
end
TestUnitRange(; start::Real, length::Integer) =
    TestUnitRange{typeof(start)}(;start=start, length=length)
Base.length(r::TestUnitRange) = r.length
Base.first(r::TestUnitRange) = r.start
Base.last(r::TestUnitRange) = @inbounds r[length(r)]
Base.show(io::IO, r::TestUnitRange) =
    print(io, "TestUnitRange{",eltype(r),"}(start=",r.start,", length=",r.length,")")
@inline function Base.getindex(r::TestUnitRange{T}, i::Int) where {T}
    @boundscheck checkbounds(Bool, r, i)
    return T(i - 1)::T + r.start
end

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

# For some numbers like `BigInt` or `BigFloat`, instances of the same value are not equal
# in the sense of `===`.
same_value_and_type(x, y) = false
same_value_and_type(x::T, y::T) where {T} = (x === y) || (x == y)

@testset "TypeUtils" begin
    @testset "Miscellaneous" begin
        # Check that TypeUtils.Unsupported cannot be instantiated.
        @test_throws Exception Unsupported()
        @test Union{Real,Unsupported} <: @inferred Unsupported(Real)
        f(x::Unsupported(Integer)) = error("this method is not yet implemented")
        @test_throws ErrorException f(8)
        f(x::Integer) = x + 1
        @test f(8) == 9
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
        @test 3         === @inferred as(Int, 3)
        @test 1         === @inferred as(Int, 1.0)
        @test Int16(1)  === @inferred as(Int16, 1.0)
        @test (1,2,3,4) === @inferred map(as(Int), (Int8(1), Int16(2), Int32(3), Int64(4)))
        @test_throws InexactError as(Int, sqrt(2))
        @test_throws Exception as(Int, π)

        @test (1,2,3)     === @inferred as(Tuple, CartesianIndex((1,2,3)))
        @test (1,2,3)     === @inferred as(NTuple{3}, CartesianIndex((1,2,3)))
        @test (2:5, -1:4) === @inferred as(Tuple, CartesianIndices((2:5, -1:4)))
        @test (2:5, -1:4) === @inferred as(NTuple{2}, CartesianIndices((2:5, -1:4)))

        @test CartesianIndex()        === @inferred as(CartesianIndex,())
        @test CartesianIndex((1,2,3)) === @inferred as(CartesianIndex, CartesianIndex((1,2,3)))
        @test CartesianIndex((1,2,3)) === @inferred as(CartesianIndex{3}, CartesianIndex((1,2,3)))
        @test CartesianIndex((1,2,3)) === @inferred as(CartesianIndex, (0x1,2,Int16(3)))
        @test CartesianIndex((1,2,3)) === @inferred as(CartesianIndex{3}, (0x1,2,Int16(3)))

        @test CartesianIndices(())           === @inferred as(CartesianIndices,())
        @test CartesianIndices((2:3,6))      === @inferred as(CartesianIndices,(2:3,6))
        @test CartesianIndices((2:3,6))      === @inferred as(CartesianIndices,CartesianIndices((2:3,6)))
        @test CartesianIndices((2:3,6,-1:4)) === @inferred as(CartesianIndices{3},(2:3,6,-1:4))
        @test CartesianIndices((2:3,6,-1:4)) === @inferred as(CartesianIndices{3},CartesianIndices((2:3,6,-1:4)))

        @test :hello === @inferred as(Symbol, "hello")
        @test "hello" == @inferred as(String, :hello)
        @test as(String, :hello) isa String
    end

    @testset "nearest()" begin
        a, b = prevfloat(1/2), nextfloat(1/2)
        for T in (BIT_INTEGERS..., BigInt)
            f = @inferred nearest(T)
            for x in (pi, 3//4, 0, 1, 4, sqrt(2), 27.0, 9f0, 0.1, true, false,
                      0.4999999999999999, 0.5, 0.5000000000000001)
                while true
                    y = round(x)
                    if T <: Bool
                        r = y > zero(y)
                    elseif !(T <: BigInt) && y <= typemin(T)
                        r = typemin(T)
                    elseif !(T <: BigInt) && y >= typemax(T)
                        r = typemax(T)
                    else
                        r = T(y)
                    end
                    @test same_value_and_type(@inferred(nearest(T,x)), @inferred(f(x)))
                    @test same_value_and_type(@inferred(nearest(T,x)), r)
                    # Next round is with -x (if it makes sense).
                    x isa Bool && break
                    x isa Irrational && break
                    x > zero(x) || break
                    x = -x
                end
            end
            @test same_value_and_type(@inferred(nearest(T, -b)), T <: Signed ? -one(T) : zero(T))
            @test same_value_and_type(@inferred(nearest(T, -a)), zero(T))
            @test same_value_and_type(@inferred(nearest(T,  a)), zero(T))
            @test same_value_and_type(@inferred(nearest(T,  b)), one(T))
            @test same_value_and_type(@inferred(nearest(T, float(zero(T)))), zero(T))
            @test same_value_and_type(@inferred(nearest(T, float(one(T)))), one(T))
            if T == BigInt
                # `typemin` and `typemax` make no sense for `BigInt`s.
                @test_throws InexactError nearest(T, -Inf)
                @test_throws InexactError nearest(T,  Inf)
            else
                @test same_value_and_type(@inferred(nearest(T, -Inf)), typemin(T))
                @test same_value_and_type(@inferred(nearest(T,  Inf)), typemax(T))
            end
            @test_throws InexactError nearest(T, NaN)
            for S in (BIT_INTEGERS..., BigInt)
                @test same_value_and_type(@inferred(nearest(T, zero(S))), zero(T))
                @test same_value_and_type(@inferred(nearest(T, one(S))), one(T))
                if T == BigInt && S != BigInt
                    @test same_value_and_type(@inferred(nearest(T, typemin(S))), T(typemin(S)))
                    @test same_value_and_type(@inferred(nearest(T, typemax(S))), T(typemax(S)))
                end
                if T != BigInt && S != BigInt
                    @test @inferred(nearest(T, typemin(S))) === (typemin(S) <= typemin(T) ? typemin(T) : T(typemin(S)))
                    @test @inferred(nearest(T, typemax(S))) === (typemax(S) >= typemax(T) ? typemax(T) : T(typemax(S)))
                end
            end
        end
        # Cases where `nearest` behaves as identity or non-integer `T`.
        @test -42    === @inferred nearest(Int, -42)
        @test -42    === @inferred nearest(Int, Int8(-42))
        @test 4.2    === @inferred nearest(Float64, 4.2)
        @test 42.0   === @inferred nearest(Float64, 42)
        @test -7.0f0 === @inferred nearest(Float32, -7)
        # Mapping.
        let vals = (1.3, -2.7, pi, 42)
            @test map(nearest(Int), vals) === map(x -> nearest(Int, x), vals)
        end
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
        @test 7 === @inferred as_array_dim(7)
        @test 9 === @inferred as_array_dim(Int16(9))
        @test 5 === @inferred as_array_dim(Base.OneTo(5))
        @test 6 === @inferred as_array_dim(-1:4)
        @test 4 === @inferred as_array_dim(Int8(0):Int8(3))
        @test 6 === @inferred as_array_dim(-1:1:4)
        @test_throws ArgumentError as_array_dim(-1:2:3) # non-unit step

        @test_throws MethodError as_array_axis(π)
        @test Base.OneTo(7) === @inferred as_array_axis(7)
        @test Base.OneTo(9) === @inferred as_array_axis(Int16(9))
        @test Base.OneTo(5) === @inferred as_array_axis(Base.OneTo(5))
        @test -1:4          === @inferred as_array_axis(-1:4)
        @test 1:4           === @inferred as_array_axis(1:4)
        @test 0:3           === @inferred as_array_axis(Int8(0):Int8(3))
        @test -1:4          === @inferred as_array_axis(-1:1:4)
        @test_throws ArgumentError as_array_axis(-1:2:3) # non-unit step

        @test_throws MethodError as_array_shape(π)
        @test_throws MethodError as_array_shape((π,))
        @test ()                                === @inferred as_array_shape()
        @test ()                                === @inferred as_array_shape(())
        @test (5,)                              === @inferred as_array_shape(5)
        @test (2:5,)                            === @inferred as_array_shape(2:5)
        @test (2,9,3)                           === @inferred as_array_shape(2,Int16(9),Int8(3))
        @test (2,9,3)                           === @inferred as_array_shape((2,Int16(9),Int8(3)))
        @test (2,3,4)                           === @inferred as_array_shape(2,3,4)
        @test (2,3,4)                           === @inferred as_array_shape(2,Base.OneTo{Int8}(3),4)
        @test (Base.OneTo(2),1:3,Base.OneTo(4)) === @inferred as_array_shape(2,Int8(1):Int8(3),4)
        @test (0:2,-4:4,-2:1)                   === @inferred as_array_shape((0:2,-4:4,-2:1))
        @test_throws ArgumentError as_array_shape((1:4, 1:2:6,)) # non-unit step

        @test_throws MethodError as_array_axes(π)
        @test_throws MethodError as_array_axes((π,))
        @test ()                                          === @inferred as_array_axes()
        @test ()                                          === @inferred as_array_axes(())
        @test (Base.OneTo(5),)                            === @inferred as_array_axes(5)
        @test (2:5,)                                      === @inferred as_array_axes(2:5)
        @test (Base.OneTo(2),Base.OneTo(9),Base.OneTo(3)) === @inferred as_array_axes(2,Int16(9),Int8(3))
        @test (Base.OneTo(2),Base.OneTo(9),Base.OneTo(3)) === @inferred as_array_axes((2,Int16(9),Int8(3)))
        @test (Base.OneTo(2),Base.OneTo(3),Base.OneTo(4)) === @inferred as_array_axes(2,3,4)
        @test (Base.OneTo(2),Base.OneTo(3),Base.OneTo(4)) === @inferred as_array_axes(2,Base.OneTo{Int8}(3),4)
        @test (Base.OneTo(2),1:3,Base.OneTo(4))           === @inferred as_array_axes(2,Int8(1):Int8(3),4)
        @test (0:2,-4:4,-2:1)                             === @inferred as_array_axes((0:2,-4:4,-2:1))
        @test_throws ArgumentError as_array_axes((1:4, 1:2:6,)) # non-unit step

        @test_throws MethodError as_array_size(π)
        @test_throws MethodError as_array_size((π,))
        @test ()      === @inferred as_array_size()
        @test ()      === @inferred as_array_size(())
        @test (5,)    === @inferred as_array_size(5)
        @test (4,)    === @inferred as_array_size(2:5)
        @test (2,9,3) === @inferred as_array_size(2,Int16(9),Int8(3))
        @test (2,9,3) === @inferred as_array_size((2,Int16(9),Int8(3)))
        @test (2,3,4) === @inferred as_array_size(2,3,4)
        @test (2,3,4) === @inferred as_array_size(2,Base.OneTo{Int8}(3),4)
        @test (2,3,4) === @inferred as_array_size(2,Int8(1):Int8(3),4)
        @test (3,9,4) === @inferred as_array_size((0:2,-4:4,-2:1))
        @test_throws ArgumentError as_array_size((1:4, 1:2:6,)) # non-unit step
    end

    @testset "new_array()" begin
        A = @inferred new_array(Int)
        @test A isa Array{Int,0}
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
        @test Int    === @inferred to_same_concrete_type(Int)
        @test UInt8  === @inferred to_same_concrete_type(UInt8, UInt8)
        @test Int8   === @inferred to_same_concrete_type(Int8, Int8, Int8)
        @test UInt16 === @inferred to_same_concrete_type(UInt8, UInt16)
        @test Int32  === @inferred to_same_concrete_type(Int8, Int16, Int32)
        @test_throws ArgumentError to_same_concrete_type()
        @test_throws ArgumentError to_same_concrete_type(AbstractUnitRange)
        @test_throws ArgumentError to_same_concrete_type(AbstractFloat, Float32)
        @test_throws ArgumentError to_same_concrete_type(String, Float32, Int16)
    end

    @testset "to_same_type()" begin
        @test ()                             === @inferred to_same_type()
        @test (pi,)                          === @inferred to_same_type(pi)
        @test ('x',)                         === @inferred to_same_type('x')
        @test (1, 2)                         === @inferred to_same_type(1, 2)
        @test (0x1, 0x2, 0x3)                === @inferred to_same_type(0x1, 0x2, 0x3)
        @test (1, 2)                         === @inferred to_same_type(0x1, 2)
        @test (Int16(1), Int16(2), Int16(3)) === @inferred to_same_type(UInt8(1), Int8(2), Int16(3))
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
        @test promote_type()                                === @inferred promote_eltype()
        @test eltype(A)                                     === @inferred promote_eltype(A)
        @test eltype(B)                                     === @inferred promote_eltype(B)
        @test eltype(C)                                     === @inferred promote_eltype(C)
        @test promote_type(eltype(A), eltype(B))            === @inferred promote_eltype(A,B)
        @test promote_type(eltype(A), eltype(B), eltype(C)) === @inferred promote_eltype(A,B,C)
    end

    @testset "convert_eltype()" begin
        # Numbers.
        let A = rand(Float64), B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa Float32
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(eltype(B), typeof(A))
        end
        # Abstract arrays.
        let A = rand(Float64, 2, 3), B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa AbstractArray{Float32,ndims(A)}
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(eltype(B), typeof(A))
            @test AbstractArray{UInt8,2} === @inferred convert_eltype(UInt8, typeof(view(A,:,1:2:3)))
        end
        # Base.OneTo
        let A = OneTo{Int32}(7), B = @inferred convert_eltype(Int16, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa OneTo{Int16}
            @test B == Int16.(A)
            @test typeof(B) === @inferred convert_eltype(eltype(B), typeof(A))
        end
        let A = OneTo(7), B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(eltype(B), typeof(A))
        end
        # UnitRange
        let A = 2:8, B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(eltype(B), typeof(A))
        end
        # TestUnitRange
        let A = TestUnitRange(;start=0.25, length=8), B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(eltype(B), typeof(A))
        end
        # OrdinalRange
        let A = 2.0:3.0:11.0, B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa AbstractRange{Float32}
            @test B == Float32.(A)
            @test typeof(B) <: @inferred convert_eltype(eltype(B), typeof(A))
        end
        # LinRange
        let A = LinRange(-2.0, 3.0, 5), B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa LinRange{Float32}
            @test B == Float32.(A)
            @test typeof(B) <: @inferred convert_eltype(eltype(B), typeof(A))
        end
        # Tuples.
        let A = (1, 2, 3) #= NTuple =#, B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa NTuple{3,Float32}
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(Float32, typeof(A))
        end
        let A = (1, 2.0, pi) #= not NTuple but Tuple =#, B = @inferred convert_eltype(Float32, A)
            @test A === @inferred convert_eltype(eltype(A), A)
            @test B isa NTuple{3,Float32}
            @test B == Float32.(A)
            @test typeof(B) === @inferred convert_eltype(Float32, typeof(A))
        end
        # Map `convert_eltype`.
        let A = rand(Float64, 2, 3), B = reshape(1:12, 3, 4)
            f = @inferred convert_eltype(Float32)
            @test isequal(f(A), convert_eltype(Float32, A))
            @test isequal(f(B), convert_eltype(Float32, B))
        end
        @test_throws ErrorException convert_eltype(Char, String)
    end

    @testset "as_eltype()" begin
        let A = rand(Float64, 3, 4, 5), B = @inferred as_eltype(Float32, A)
            @test A === @inferred as_eltype(eltype(A), A)
            @test eltype(B) === Float32
            @test length(B) === length(A)
            @test size(B) === size(A)
            @test axes(B) === axes(A)
            @test IndexStyle(B) === IndexStyle(A)
            @test parent(B) === A
            @test B == Float32.(A)
            A[1,2,3] = -7
            @test B[1,2,3] === Float32(-7)
            B[1,2,3] = 19
            @test A[1,2,3] == 19
            C = copyto!(@inferred(similar(B)), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), size(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), axes(B))), B)
            @test C == B
        end
        let A = view(rand(Float64, 3, 4, 5), :, 2, :), B = @inferred as_eltype(Float32, A)
            @test A === @inferred as_eltype(eltype(A), A)
            @test eltype(B) == Float32
            @test length(B) === length(A)
            @test size(B) === size(A)
            @test axes(B) === axes(A)
            @test IndexStyle(B) === IndexStyle(A)
            @test parent(B) === A
            @test B == Float32.(A)
            A[2,3] = -7
            @test B[2,3] === Float32(-7)
            B[2,3] = 19
            @test A[2,3] == 19
            C = copyto!(@inferred(similar(B)), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), size(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), axes(B))), B)
            @test C == B
        end
        let A = 1:5, B = @inferred as_eltype(Float32, A)
            @test A === @inferred as_eltype(eltype(A), A)
            @test length(B) === length(A)
            @test size(B) === size(A)
            @test axes(B) === axes(A)
            @test IndexStyle(B) === IndexStyle(A)
            @test parent(B) === A
            @test B == Float32.(A)
            C = copyto!(@inferred(similar(B)), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), size(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), axes(B))), B)
            @test C == B
        end
        let A = 2.0:3.0:11.0, B = @inferred as_eltype(Float32, A)
            @test A === @inferred as_eltype(eltype(A), A)
            @test length(B) === length(A)
            @test size(B) === size(A)
            @test axes(B) === axes(A)
            @test IndexStyle(B) === IndexStyle(A)
            @test parent(B) === A
            @test B == Float32.(A)
            C = copyto!(@inferred(similar(B)), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), size(B))), B)
            @test C == B
            C = copyto!(@inferred(similar(B, eltype(B), axes(B))), B)
            @test C == B
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
        @test Base.return_types(f) == [return_type(f)]
        @test Base.return_types(g, Tuple) == [return_type(g)]
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
        @test BareNumber === @inferred bare_type()

        # bare_type for values
        @test Float64      === @inferred bare_type(1.0)
        @test Float32      === @inferred bare_type(Float32)
        @test Complex{Int} === @inferred bare_type(Complex(2,3))
        @test typeof(NaN)  === @inferred bare_type(NaN)
        @test typeof(π)    === @inferred bare_type(π)
        @test typeof(3//4) === @inferred bare_type(3//4)
        @test_throws ErrorException bare_type("hello")

        # bare_type for types
        @test Real         === @inferred bare_type(Real)
        @test Integer      === @inferred bare_type(Integer)
        @test Float32      === @inferred bare_type(Float32)
        @test BigFloat     === @inferred bare_type(BigFloat)
        @test Complex{Int} === @inferred bare_type(Complex{Int})
        @test typeof(π)    === @inferred bare_type(typeof(π))
        @test typeof(3//4) === @inferred bare_type(typeof(3//4))
        @test_throws ErrorException bare_type(AbstractString)

        # bare_type with multiple arguments
        @test Float32                                               === @inferred bare_type(1, 0f0)
        @test promote_type(Int, typeof(pi))                         === @inferred bare_type(Int, pi)
        @test promote_type(Int, typeof(pi), Float64)                === @inferred bare_type(4, pi, 1.0)
        @test promote_type(Int, Int8, Float32)                      === @inferred bare_type(Int, Int8, Float32)
        @test promote_type(Int, Int8, Int16, Float32)               === @inferred bare_type(Int, Int8, Int16, Float32)
        @test promote_type(Int, Int8, Int16, Float32, Bool)         === @inferred bare_type(Int, Int8, Int16, Float32, Bool)
        @test promote_type(Int, Int8, Int16, Float32, Bool, UInt16) === @inferred bare_type(Int, Int8, Int16, Float32, Bool, UInt16)

        # default implementation
        @test Float32 === @inferred bare_type(MyNumber(1.2f0))
        @test Int16   === @inferred bare_type(MyNumber{Int16})

        # convert_bare_type
        @test -1        === @inferred convert_bare_type(Int, -1)
        @test 2         === @inferred convert_bare_type(Int, 2.0)
        @test 2.0f0     === @inferred convert_bare_type(Float32, 2.0)
        @test Int16(12) === @inferred convert_bare_type(MyNumber{Int16}, 12.0)
        @test_throws ErrorException convert_bare_type(Int, "oups!")
        for x in (missing, nothing, undef)
            @test x         === @inferred convert_bare_type(Int8, x)
            @test typeof(x) === @inferred convert_bare_type(Int8, typeof(x))
        end
        @test Int16                    === @inferred convert_bare_type(Int16, Int32)
        @test Int16                    === @inferred convert_bare_type(Int16, Complex{Int32})
        @test Complex{Int16}           === @inferred convert_bare_type(Complex{Int16}, Int32)
        @test Complex{Int16}           === @inferred convert_bare_type(Complex{Int16}, Complex{Int32})
        @test MyNumber{Int16}          === @inferred convert_bare_type(MyNumber{Int16}, MyNumber{Int32})
        @test MyNumber{Int16}          === @inferred convert_bare_type(Int16, MyNumber{Complex{Int32}})
        @test MyNumber{Complex{Int16}} === @inferred convert_bare_type(Complex{Int16}, MyNumber{Int32})
        @test_throws ErrorException convert_bare_type(Int, String)
        let vals = (1.3, -2.7, pi, 42, u"35GHz")
            @test map(x -> convert_bare_type(Float32, x), vals) === @inferred map(convert_bare_type(Float32), vals)
        end

        # real_type with no argument
        @test Real === @inferred real_type()

        # real_type for values
        @test Float64      === @inferred real_type(1.0)
        @test Float32      === @inferred real_type(Float32)
        @test Int          === @inferred real_type(Complex(2,3))
        @test typeof(NaN)  === @inferred real_type(NaN)
        @test typeof(π)    === @inferred real_type(π)
        @test typeof(3//4) === @inferred real_type(3//4)
        @test_throws ErrorException real_type("hello")

        # real_type for types
        @test Real         === @inferred real_type(Real)
        @test Integer      === @inferred real_type(Integer)
        @test Float32      === @inferred real_type(Float32)
        @test BigFloat     === @inferred real_type(BigFloat)
        @test Int          === @inferred real_type(Complex{Int})
        @test typeof(π)    === @inferred real_type(typeof(π))
        @test typeof(3//4) === @inferred real_type(typeof(3//4))
        @test_throws ErrorException real_type(AbstractString)

        # real_type with multiple arguments
        @test Float32                                               === @inferred real_type(1, 0f0)
        @test promote_type(Int, typeof(pi))                         === @inferred real_type(Int, pi)
        @test promote_type(Int, typeof(pi), Float64)                === @inferred real_type(4, pi, 1.0)
        @test promote_type(Int, Int8, Float32)                      === @inferred real_type(Int, Int8, Float32)
        @test promote_type(Int, Int8, Int16, Float32)               === @inferred real_type(Int, Int8, Int16, Float32)
        @test promote_type(Int, Int8, Int16, Float32, Bool)         === @inferred real_type(Int, Int8, Int16, Float32, Bool)
        @test promote_type(Int, Int8, Int16, Float32, Bool, UInt16) === @inferred real_type(Int, Int8, Int16, Float32, Bool, UInt16)
        @test promote_type(Int, Int8, Int16, Float32, Bool, UInt16) === @inferred real_type(Int, Int8, Int16, Complex{Float32}, Bool, UInt16)

        # default implementation
        @test Float32 === @inferred real_type(MyNumber(1.2f0))
        @test Int16   === @inferred real_type(MyNumber{Complex{Int16}})

        # convert_real_type
        @test -1    === @inferred convert_real_type(Int, -1)
        @test 2     === @inferred convert_real_type(Int, 2.0)
        @test 2.0f0 === @inferred convert_real_type(Complex{Float32}, 2.0)
        let x = 2.1f0
            @test x === @inferred convert_real_type(typeof(x), x)
        end
        let x = 2 + 3im
            @test x === @inferred convert_real_type(typeof(x), x)
            @test x === @inferred convert_real_type(real(typeof(x)), x)
        end
        @test Complex{Int16}(2, 3)    === @inferred convert_real_type(Complex{Int16}, 2 + 3im)
        @test Complex{Float32}(2, -1) === @inferred convert_real_type(Float32, 2.0 - 1.0im)
        @test Int16(12)               === @inferred convert_real_type(MyNumber{Int16}, 12.0)
        @test_throws ErrorException convert_real_type(Int, "oups!")
        for x in (missing, nothing, undef)
            @test x         === @inferred convert_real_type(Int8, x)
            @test typeof(x) === @inferred convert_real_type(Int8, typeof(x))
        end
        @test Int16                    === @inferred convert_real_type(Int16, Int32)
        @test Complex{Int16}           === @inferred convert_real_type(Int16, Complex{Int32})
        @test Int16                    === @inferred convert_real_type(Complex{Int16}, Int32)
        @test Complex{Int16}           === @inferred convert_real_type(Complex{Int16}, Complex{Int32})
        @test MyNumber{Int16}          === @inferred convert_real_type(MyNumber{Int16}, MyNumber{Int32})
        @test MyNumber{Complex{Int16}} === @inferred convert_real_type(Int16, MyNumber{Complex{Int32}})
        @test MyNumber{Int16}          === @inferred convert_real_type(Complex{Int16}, MyNumber{Int32})
        @test_throws ErrorException convert_real_type(Int, String)
        let vals = (1.3, -2.7, pi, Complex(-3,42), u"35GHz")
            @test map(x -> convert_real_type(Float32, x), vals) === @inferred map(convert_real_type(Float32), vals)
        end

        # floating_point_type
        @test AbstractFloat === @inferred floating_point_type()
        @test float(Int16)  === @inferred floating_point_type(Int16)
        @test Float32       === @inferred floating_point_type(Int16, Float32)
        @test Float64       === @inferred floating_point_type(Int16, Float32, Complex{Float64})
        @test Float32       === @inferred floating_point_type(Int16, Float32, Complex{Float16}, 0x1)
        @test Float32       === @inferred floating_point_type(Int16, Float32, Complex{Float16}, 0x1, -42)

        # convert_floating_point_type
        @test -1.0  === @inferred convert_floating_point_type(Int, -1)
        @test 2.0   === @inferred convert_floating_point_type(Int, 2.0)
        @test 2.0f0 === @inferred convert_floating_point_type(Complex{Float32}, 2.0)
        let x = 2.1f0
            @test x === @inferred convert_floating_point_type(typeof(x), x)
        end
        let x = 2 + 3im
            @test float(x) === @inferred convert_floating_point_type(typeof(x), x)
            @test float(x) === @inferred convert_floating_point_type(real(typeof(x)), x)
        end
        @test Complex{Float64}(2, 3)  === @inferred convert_floating_point_type(Complex{Int16}, 2 + 3im)
        @test Complex{Float32}(2, -1) === @inferred convert_floating_point_type(Float32, 2.0 - 1.0im)
        @test Float64(12)             === @inferred convert_floating_point_type(MyNumber{Int16}, 12.0)
        @test_throws ErrorException convert_floating_point_type(Int, "oups!")
        for x in (missing, nothing, undef)
            @test x         === @inferred convert_floating_point_type(Float32, x)
            @test typeof(x) === @inferred convert_floating_point_type(Float32, typeof(x))
        end
        @test Float32                    === @inferred convert_floating_point_type(Float32, Float64)
        @test Complex{Float32}           === @inferred convert_floating_point_type(Float32, Complex{Float64})
        @test Float32                    === @inferred convert_floating_point_type(Complex{Float32}, Float64)
        @test Complex{Float32}           === @inferred convert_floating_point_type(Complex{Float32}, Complex{Float64})
        @test MyNumber{Float32}          === @inferred convert_floating_point_type(MyNumber{Float32}, MyNumber{Float64})
        @test MyNumber{Complex{Float32}} === @inferred convert_floating_point_type(Float32, MyNumber{Complex{Float64}})
        @test MyNumber{Float32}          === @inferred convert_floating_point_type(Complex{Float32}, MyNumber{Float64})
        @test_throws ErrorException convert_floating_point_type(Int, String)
        let vals = (1.3, -2.7, pi, Complex(-3,42), u"35GHz")
            @test map(x -> convert_floating_point_type(Float32, x), vals) === @inferred map(convert_floating_point_type(Float32), vals)
        end

        # assert_floating_point
        @test !assert_floating_point(Bool, 1)
        @test !assert_floating_point(Bool, π)
        @test !assert_floating_point(Bool, 1//2)
        @test !assert_floating_point(Bool, "hello")
        @test !assert_floating_point(Bool, [0x1, 0x2])
        @test !assert_floating_point(Bool, (0x1, -1))
        @test  assert_floating_point(Bool, 1.0)
        @test  assert_floating_point(Bool, complex(-2f0,3f0))
        @test  assert_floating_point(Bool, [0.0f0, 2.0f0])
        @test !assert_floating_point(Bool, Int)
        @test !assert_floating_point(Bool, AbstractFloat)
        @test  assert_floating_point(Bool, Float16)
        @test  assert_floating_point(Bool, Float32)
        @test  assert_floating_point(Bool, Float64)
        @test  assert_floating_point(Bool, Complex{Float32})
        @test  assert_floating_point(Bool, AbstractRange{Float64})
        # only ok for homogeneous tuples
        @test !assert_floating_point(Bool, (-1.0, 2.0f0))
        @test  assert_floating_point(Bool, (-1.0, 2.0))
        @test !assert_floating_point(Bool, Tuple{Float64,Float32})
        @test  assert_floating_point(Bool, Tuple{Float32,Float32})
        #
        x, y = 2, 3.0
        @test assert_floating_point(y) === nothing
        @test assert_floating_point("y", y) === nothing
        @test_throws Exception assert_floating_point(x)
        @test_throws Exception assert_floating_point(:x, x)
        @test_throws Exception @assert_floating_point x y

        # unitless
        @test bare_type(Real)         === @inferred unitless(Real)
        @test bare_type(Integer)      === @inferred unitless(Integer)
        @test bare_type(Float32)      === @inferred unitless(Float32)
        @test bare_type(BigFloat)     === @inferred unitless(BigFloat)
        @test bare_type(Complex{Int}) === @inferred unitless(Complex{Int})
        @test bare_type(typeof(3//4)) === @inferred unitless(typeof(3//4))
        @test bare_type(typeof(π))    === @inferred unitless(typeof(π))
        @test 17.0                    === @inferred unitless(17.0)
        @test 17.0f0                  === @inferred unitless(17.0f0)
        @test 17                      === @inferred unitless(17)
        @test Int16(17)               === @inferred unitless(Int16(17))
        @test true                    === @inferred unitless(true)
        @test false                   === @inferred unitless(false)
        @test 3//4                    === @inferred unitless(3//4)
        @test π                       === @inferred unitless(π)
        @test 1.2f0                   === @inferred unitless(MyNumber(1.2f0))
        @test Int16                   === @inferred unitless(MyNumber{Int16})

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

    @testset "get_precision" begin

        @test @inferred(get_precision(Symbol                  )) == AbstractFloat
        @test @inferred(get_precision(String                  )) == AbstractFloat
        @test @inferred(get_precision(Bool                    )) == AbstractFloat
        @test @inferred(get_precision(UInt8                   )) == AbstractFloat
        @test @inferred(get_precision(UInt16                  )) == AbstractFloat
        @test @inferred(get_precision(UInt32                  )) == AbstractFloat
        @test @inferred(get_precision(UInt64                  )) == AbstractFloat
        @test @inferred(get_precision(UInt128                 )) == AbstractFloat
        @test @inferred(get_precision(Int8                    )) == AbstractFloat
        @test @inferred(get_precision(Int16                   )) == AbstractFloat
        @test @inferred(get_precision(Int32                   )) == AbstractFloat
        @test @inferred(get_precision(Int64                   )) == AbstractFloat
        @test @inferred(get_precision(Int128                  )) == AbstractFloat
        @test @inferred(get_precision(BigInt                  )) == AbstractFloat
        @test @inferred(get_precision(Float16                 )) == Float16
        @test @inferred(get_precision(Float32                 )) == Float32
        @test @inferred(get_precision(Float64                 )) == Float64
        @test @inferred(get_precision(BigFloat                )) == BigFloat
        @test @inferred(get_precision(Complex{Int}            )) == AbstractFloat
        @test @inferred(get_precision(Complex{Float16}        )) == Float16
        @test @inferred(get_precision(Complex{Float32}        )) == Float32
        @test @inferred(get_precision(Complex{Float64}        )) == Float64
        @test @inferred(get_precision(Complex{BigFloat}       )) == BigFloat
        @test @inferred(get_precision(Array{Int,3}            )) == AbstractFloat
        @test @inferred(get_precision(Array{Float64}          )) == Float64
        @test @inferred(get_precision(Vector{BigFloat}        )) == BigFloat
        @test @inferred(get_precision(Matrix{Complex{Float16}})) == Float16
        @test @inferred(get_precision(UnitRange{BigFloat}     )) == BigFloat

        @test @inferred(get_precision(:symbol                  )) == AbstractFloat
        @test @inferred(get_precision("string"                 )) == AbstractFloat
        @test @inferred(get_precision(true                     )) == AbstractFloat
        @test @inferred(get_precision(false                    )) == AbstractFloat
        @test @inferred(get_precision(3                        )) == AbstractFloat
        @test @inferred(get_precision(0x03                     )) == AbstractFloat
        @test @inferred(get_precision(3                        )) == AbstractFloat
        @test @inferred(get_precision(big(3)                   )) == AbstractFloat
        @test @inferred(get_precision(3//2                     )) == AbstractFloat
        @test @inferred(get_precision(π                        )) == AbstractFloat
        @test @inferred(get_precision(Float16(3)               )) == Float16
        @test @inferred(get_precision(3.0f0                    )) == Float32
        @test @inferred(get_precision(3.0                      )) == Float64
        @test @inferred(get_precision(big(3.0)                 )) == BigFloat
        @test @inferred(get_precision(1 + 2im                  )) == AbstractFloat
        @test @inferred(get_precision(1.0f0 + 2.0f0im          )) == Float32
        @test @inferred(get_precision(1.0 + 2.0im              )) == Float64
        @test @inferred(get_precision(-2:5                     )) == AbstractFloat
        @test @inferred(get_precision(-2.0f0:5.0f0             )) == Float32
        @test @inferred(get_precision(-2.0:5.0                 )) == Float64
        @test @inferred(get_precision((true, ("", 2), :sym)    )) == AbstractFloat
        @test @inferred(get_precision(("", 2, big(1.0))        )) == BigFloat
        @test @inferred(get_precision((; a=false, b="", c=-2.1))) == Float64
        @test @inferred(get_precision((; a=false, b="", c=π)   )) == AbstractFloat

        A = ones(Bool, 2,3,4)
        @test @inferred(get_precision(A)) == AbstractFloat
        @test @inferred(get_precision(typeof(A))) == AbstractFloat

        A = ComplexF32.([9+1im 2-3im 1; 0 7 1; 0 0 4])
        @test @inferred(get_precision(A)) == Float32
        @test @inferred(get_precision(typeof(A))) == Float32
        B = adjoint(A)
        @test @inferred(get_precision(B)) == Float32
        @test @inferred(get_precision(typeof(B))) == Float32
        B = Diagonal(A)
        @test @inferred(get_precision(B)) == Float32
        @test @inferred(get_precision(typeof(B))) == Float32
        B = Hermitian(A)
        @test @inferred(get_precision(B)) == Float32
        @test @inferred(get_precision(typeof(B))) == Float32
        B = qr(A)
        @test @inferred(get_precision(B)) == Float32
        @test @inferred(get_precision(typeof(B))) == Float32

    end

    @testset "adapt_precision" begin
        # Precision converter.
        @test_throws Exception adapt_precision(String, (1,2,:x))
        @test_throws Exception adapt_precision(Symbol)
        @test TypeUtils.default_precision == Float64
        f = @inferred adapt_precision(AbstractFloat)
        @test f === @inferred adapt_precision(TypeUtils.default_precision)
        @test 7.0 === @inferred f(7)
        @test 7.0 === @inferred f(7.0f0)
        @test 7.0 === @inferred f(big(7))
        @test 7.0 === @inferred f(big(7.0))
        @test -3.0u"kg/s" === @inferred f(-3u"kg/s")

        @testset "adapt_precision($T, x)" for T in (AbstractFloat, Float16, Float32, Float64, BigFloat)
            if isconcretetype(T)
                @test T <: TypeUtils.Precision
            end

            S = isconcretetype(T) ? T : TypeUtils.default_precision

            @test @inferred(adapt_precision(T, Symbol    )) === Symbol
            @test @inferred(adapt_precision(T, String    )) === String
            @test @inferred(adapt_precision(T, Bool      )) === S
            @test @inferred(adapt_precision(T, UInt8     )) === S
            @test @inferred(adapt_precision(T, UInt16    )) === S
            @test @inferred(adapt_precision(T, UInt32    )) === S
            @test @inferred(adapt_precision(T, UInt64    )) === S
            @test @inferred(adapt_precision(T, UInt128   )) === S
            @test @inferred(adapt_precision(T, Int8      )) === S
            @test @inferred(adapt_precision(T, Int16     )) === S
            @test @inferred(adapt_precision(T, Int32     )) === S
            @test @inferred(adapt_precision(T, Int64     )) === S
            @test @inferred(adapt_precision(T, Int128    )) === S
            @test @inferred(adapt_precision(T, BigInt    )) === S
            @test @inferred(adapt_precision(T, Float16   )) === S
            @test @inferred(adapt_precision(T, Float32   )) === S
            @test @inferred(adapt_precision(T, Float64   )) === S

            str = "string"
            @test @inferred(adapt_precision(T, :symbol   )) === :symbol
            @test @inferred(adapt_precision(T, str       )) === str # same object
            @test @inferred(adapt_precision(T, true      )) ≗ one(S)
            @test @inferred(adapt_precision(T, false     )) ≗ zero(S)
            @test @inferred(adapt_precision(T, 0x03      )) ≗ S(3)
            @test @inferred(adapt_precision(T, 3         )) ≗ S(3)
            @test @inferred(adapt_precision(T, big(3)    )) ≗ S(3)
            @test @inferred(adapt_precision(T, 3//2      )) ≗ S(3//2)
            @test @inferred(adapt_precision(T, π         )) ≗ S(π)
            @test @inferred(adapt_precision(T, Float16(3))) ≗ S(3)
            @test @inferred(adapt_precision(T, 3.0f0     )) ≗ S(3)
            @test @inferred(adapt_precision(T, 3.0       )) ≗ S(3)
            @test @inferred(adapt_precision(T, big(3.0)  )) ≗ S(3)

            # Tests with an array.
            A = reshape(-3:20, 2,3,4)
            B = @inferred adapt_precision(T, A)
            @test eltype(B) === S
            @test axes(B) == axes(A)
            @test B == A
            @test typeof(B) <: @inferred adapt_precision(T, typeof(A))
            @test typeof(B) === @inferred adapt_precision(T, typeof(B))
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            # Tests with a real rectangular matrix.
            A = [4 -1 0; 2 3 -1]
            AAt = A*A'
            @test S == @inferred adapt_precision(T, eltype(A))
            B = @inferred adapt_precision(T, A)
            @test eltype(B) == S
            @test axes(B) == axes(A)
            @test B == A
            @test typeof(B) <: @inferred adapt_precision(T, typeof(A))
            @test typeof(B) == @inferred adapt_precision(T, typeof(B))
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = adjoint(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Adjoint
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = transpose(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Transpose
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = Symmetric(AAt)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Symmetric
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = Diagonal(AAt)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Diagonal
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = Bidiagonal(AAt, :U)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Bidiagonal
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = Tridiagonal(AAt)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Tridiagonal
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = cholesky(AAt)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Cholesky
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B ≃ C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = qr(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Factorization
                @test eltype(B) == S
                @test axes(B) == axes(C)
                if T != Float16 && T != BigFloat
                    @test B ≃ C
                end
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = svd(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Factorization
                @test eltype(B) == S
                @test axes(B) == axes(C)
                @test B ≃ C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            # Tests with a complex matrix.
            A = [9+1im 2-3im 1; 0 7 1; 0 0 4] # eltype(A) = Complex{Int}
            @test Complex{S} == @inferred adapt_precision(T, eltype(A))
            B = @inferred adapt_precision(T, A)
            @test eltype(B) == Complex{S}
            @test axes(B) == axes(A)
            @test B == A
            @test typeof(B) <: @inferred adapt_precision(T, typeof(A))
            @test typeof(B) == @inferred adapt_precision(T, typeof(B))
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = adjoint(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Adjoint
                @test eltype(B) == Complex{S}
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = transpose(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Transpose
                @test eltype(B) == Complex{S}
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = Diagonal(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Diagonal
                @test eltype(B) == Complex{S}
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = Hermitian(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Hermitian
                @test eltype(B) == Complex{S}
                @test axes(B) == axes(C)
                @test B == C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = qr(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Factorization
                @test eltype(B) == Complex{S}
                @test axes(B) == axes(C)
                if T != Float16 && T != BigFloat
                    @test B ≃ C
                end
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object

            C = svd(A)
            B = @inferred adapt_precision(T, C)
            if real(eltype(C)) == T
                @test B === C # must be same object
            else
                @test typeof(B) <: Factorization
                @test eltype(B) == Complex{S}
                @test axes(B) == axes(C)
                @test B ≃ C
            end
            @test B === @inferred adapt_precision(T, B) # same precision yields identical object
        end

        @testset "adapt_precision(T, x::Tuple)" begin
            @test adapt_precision(
                Float32, (
                    true,  1,     0x07,  ("hello", 1.0,   3.0   - 2.0im,           π))) === (
                        1.0f0, 1.0f0, 7.0f0, ("hello", 1.0f0, 3.0f0 - 2.0f0im, Float32(π)))
            x = (; a = false, b = :symbol, d = [-1,   2])
            y = (; a = 0.0,   b = :symbol, d = [-1.0, 2.0])
            z = @inferred adapt_precision(Float64, x)
            @test z == y
            @test typeof(z.d) == typeof(y.d)
            # a very long tuple
            x = ([(isodd(x) ? Int8(   x) : string(x)) for x in -20:30]...,)
            y = ([(isodd(x) ? Float32(x) : string(x)) for x in -20:30]...,)
            z = adapt_precision(Float32, x)
            @test z == y
            @test get_precision(x) === AbstractFloat
            @test get_precision(y) === get_precision(z) === Float32
        end
    end

    @testset "LinearAlgebra" begin
        A = ComplexF32.([9+1im 2-3im 1; 0 7 1; 0 0 4])
        AAt = A*A'
        T = ComplexF64 # a wider type than eltype(A) for exact conversion

        B = transpose(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Transpose{T}
        @test C == B

        B = adjoint(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Adjoint{T}
        @test C == B

        B = Diagonal(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Diagonal{T}
        @test C == B

        B = Bidiagonal(A, :U)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Bidiagonal{T}
        @test C == B

        B = Bidiagonal(A, :L)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Bidiagonal{T}
        @test C == B

        B = Tridiagonal(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Tridiagonal{T}
        @test C == B

        B = Hermitian(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Hermitian{T}
        @test C == B

        B = Symmetric(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Symmetric{T}
        @test C == B

        B = cholesky(AAt)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Cholesky{T}
        @test C.factors == B.factors

        pivot = VERSION < v"1.8.0-beta1" ? Val(true) : RowMaximum()
        B = cholesky(AAt, pivot)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: CholeskyPivoted{T}
        @test C.factors == B.factors

        B = svd(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: SVD{T}
        @test C.U == B.U
        @test C.S == B.S
        @test C.Vt == B.Vt

        B = qr(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: Factorization{T}
        @test C.Q == B.Q
        @test C.R == B.R

        pivot = isdefined(LinearAlgebra, :PivotingStrategy) ? ColumnNorm() : Val(true)
        B = qr(A, pivot)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: QRPivoted{T}
        @test C.P == B.P
        @test C.Q == B.Q
        @test C.R == B.R
        @test C.p == B.p

        pivot = isdefined(LinearAlgebra, :PivotingStrategy) ? NoPivot() : Val(false)
        B = lu(A, pivot)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: LU{T}
        @test C.L == B.L
        @test C.U == B.U
        @test C.P == B.P
        @test C.p == B.p

        pivot = isdefined(LinearAlgebra, :PivotingStrategy) ? RowMaximum() : Val(true)
        B = lu(A, pivot)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        C = @inferred convert_eltype(T, B)
        @test typeof(C) <: LU{T}
        @test C.L == B.L
        @test C.U == B.U
        @test C.P == B.P
        @test C.p == B.p

        B = hessenberg(A)
        @test eltype(B) == eltype(A) != T
        @test @inferred(convert_eltype(eltype(B), B)) === B
        @test_throws Exception convert_eltype(T, B)
    end

    @testset "Unitful quantities" begin
        # bare_type for values
        @test Float64 === @inferred bare_type(u"2.0m/s")
        @test Int     === @inferred bare_type(u"35GHz")

        # bare_type for types
        @test Float64 === @inferred bare_type(typeof(u"2.0m/s"))
        @test Int     === @inferred bare_type(typeof(u"35GHz"))

        # bare_type with multiple arguments
        @test Float64          === @inferred bare_type(u"2.0m/s", u"35GHz")
        @test Float32          === @inferred bare_type(1, u"2.0f0m/s", u"35GHz")
        @test Complex{Float32} === @inferred bare_type(1, u"2.0f0m/s", u"35GHz", Complex{Int8}(11))

        # convert_bare_type
        @test u"2.0m/s"    === @inferred convert_bare_type(Float64, u"2.0m/s")
        @test u"2m/s"      === @inferred convert_bare_type(Int,     u"2.0m/s")
        @test u"35.0f0GHz" === @inferred convert_bare_type(Float32, u"35GHz")
        let T = typeof(u"3.5GHz")
            @test T === @inferred convert_bare_type(Float64, T)
            for x in (missing, nothing, undef)
                @test x         === @inferred convert_bare_type(T, x)
                @test typeof(x) === @inferred convert_bare_type(T, typeof(x))
            end
        end
        let u = u"km/s"
            @test typeof(one(Int16)*u)          === @inferred convert_bare_type(Int16, typeof(one(Int32)*u))
            @test typeof(one(Int16)*u)          === @inferred convert_bare_type(Int16, typeof(one(Complex{Int32})*u))
            @test typeof(one(Complex{Int16})*u) === @inferred convert_bare_type(Complex{Int16}, typeof(one(Int32)*u))
            @test typeof(one(Complex{Int16})*u) === @inferred convert_bare_type(Complex{Int16}, typeof(one(Complex{Int32})*u))
        end

        # real_type for values
        @test Float64 === @inferred real_type(u"2.0m/s")
        @test Int     === @inferred real_type(u"35GHz")

        # real_type for types
        @test Float64 === @inferred real_type(typeof(u"2.0m/s"))
        @test Int     === @inferred real_type(typeof(u"35GHz"))

        # real_type with multiple arguments
        @test Float64 === @inferred real_type(u"2.0m/s", u"35GHz")
        @test Float32 === @inferred real_type(1, u"2.0f0m/s", u"35GHz")
        @test Float32 === @inferred real_type(1, u"2.0f0m/s", u"35GHz", Complex{Int8}(11))

        # convert_real_type
        @test u"2.0m/s"    === @inferred convert_real_type(Float64, u"2.0m/s")
        @test u"2m/s"      === @inferred convert_real_type(Int,     u"2.0m/s")
        @test u"35.0f0GHz" === @inferred convert_real_type(Float32, u"35GHz")
        let T = typeof(u"3.5GHz")
            @test T === @inferred convert_real_type(Float64, T)
            for x in (missing, nothing, undef)
                @test x         === @inferred convert_real_type(T, x)
                @test typeof(x) === @inferred convert_real_type(T, typeof(x))
            end
        end
        let u = u"km/s"
            @test typeof(one(Int16)*u)          === @inferred convert_real_type(Int16, typeof(one(Int32)*u))
            @test typeof(one(Complex{Int16})*u) === @inferred convert_real_type(Int16, typeof(one(Complex{Int32})*u))
            @test typeof(one(Int16)*u)          === @inferred convert_real_type(Complex{Int16}, typeof(one(Int32)*u))
            @test typeof(one(Complex{Int16})*u) === @inferred convert_real_type(Complex{Int16}, typeof(one(Complex{Int32})*u))
        end

        # convert_floating_point_type
        @test u"2.0m/s"    === @inferred convert_floating_point_type(Float64, u"2.0m/s")
        @test u"2.0m/s"    === @inferred convert_floating_point_type(Int, u"2.0m/s")
        @test u"35.0f0GHz" === @inferred convert_floating_point_type(Float32, u"35GHz")
        let T = typeof(u"3.5GHz")
            for x in (missing, nothing, undef)
                @test x         === @inferred convert_floating_point_type(T, x)
                @test typeof(x) === @inferred convert_floating_point_type(T, typeof(x))
            end
        end
        let u = u"km/s"
            @test typeof(float(one(Int16))*u)          === @inferred convert_floating_point_type(Int16, typeof(one(Int32)*u))
            @test typeof(float(one(Complex{Int16}))*u) === @inferred convert_floating_point_type(Int16, typeof(one(Complex{Int32})*u))
            @test typeof(float(one(Int16))*u)          === @inferred convert_floating_point_type(Complex{Int16}, typeof(one(Int32)*u))
            @test typeof(float(one(Complex{Int16}))*u) === @inferred convert_floating_point_type(Complex{Int16}, typeof(one(Complex{Int32})*u))
        end

        # get_precision
        @test AbstractFloat === @inferred get_precision(typeof(         35  * u"GHz"))
        @test Float16       === @inferred get_precision(typeof(Float16( 35) * u"GHz"))
        @test Float32       === @inferred get_precision(typeof(Float32( 35) * u"GHz"))
        @test Float64       === @inferred get_precision(typeof(Float64( 35) * u"GHz"))
        @test AbstractFloat === @inferred get_precision(typeof(BigInt(  35) * u"GHz"))
        @test BigFloat      === @inferred get_precision(typeof(BigFloat(35) * u"GHz"))
        #
        @test AbstractFloat === @inferred get_precision(                35  * u"GHz" )
        @test Float16       === @inferred get_precision(       Float16( 35) * u"GHz" )
        @test Float32       === @inferred get_precision(       Float32( 35) * u"GHz" )
        @test Float64       === @inferred get_precision(       Float64( 35) * u"GHz" )
        @test AbstractFloat === @inferred get_precision(       BigInt(  35) * u"GHz" )
        @test BigFloat      === @inferred get_precision(       BigFloat(35) * u"GHz" )

        @testset "adapt_precision($T, x::Quantity)" for T in (AbstractFloat, Float16, Float32, Float64, BigFloat)
            S = isconcretetype(T) ? T : TypeUtils.default_precision
            u = u"GHz"
            v = -35
            r = S(v)*u
            @test typeof(r) === @inferred adapt_precision(T, typeof(         v  * u))
            @test typeof(r) === @inferred adapt_precision(T, typeof(Float16( v) * u))
            @test typeof(r) === @inferred adapt_precision(T, typeof(Float32( v) * u))
            @test typeof(r) === @inferred adapt_precision(T, typeof(Float64( v) * u))
            @test typeof(r) === @inferred adapt_precision(T, typeof(BigInt(  v) * u))
            @test typeof(r) === @inferred adapt_precision(T, typeof(BigFloat(v) * u))
            #
            @test r ≗ @inferred adapt_precision(T,          v  * u)
            @test r ≗ @inferred adapt_precision(T, Float16( v) * u)
            @test r ≗ @inferred adapt_precision(T, Float32( v) * u)
            @test r ≗ @inferred adapt_precision(T, Float64( v) * u)
            @test r ≗ @inferred adapt_precision(T, BigInt(  v) * u)
            @test r ≗ @inferred adapt_precision(T, BigFloat(v) * u)
        end

        # unitless
        @test 17      === @inferred unitless(u"17GHz")
        @test Float32 === @inferred unitless(typeof(u"2.0f0m/s"))

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

        @testset "`is_signed(x::$(typeof(x))`" for (x,r) in (
            (true, false), (false, false),
            (Int8(0), true), (UInt8(1), false),
            (Int16(0), true), (UInt16(1), false),
            (Int32(0), true), (UInt32(1), false),
            (Int64(0), true), (UInt64(1), false),
            (0.0, true), (1.0f0, true),
            (big(3), true), (big(1.0), true),
            (false//true, false), (0x2//0x3, false), (2//3, true),
            (complex(true,false), false), (complex(1//2,2//3), true),
            (complex(0x01//0x02,0x02//0x03), false), (1.0 - 4.0im, true),
            (1.3u"m/s", true), (0x0007*u"kg", false),
            )
            @test is_signed(typeof(x)) == r
            @test is_signed(x) === r
        end
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
        @test vals === @inferred destructure(Tuple, obj)
        vec = destructure(Vector, obj)
        @test vec == collect(vals)
        @test (vec...,) === vals
        vec = destructure(Vector{Float32}, obj)
        @test eltype(vec) == Float32
        @test length(vec) == length(vals)
        @test vec ≈ collect(vals)
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
                # FIXME Not all versions of Julia correctly infer the result.
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
                # FIXME Not all versions of Julia correctly infer the result.
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
    if VERSION ≥ v"1.6"
        include("aqua.jl")
    end
end

end # module TestingTypeUtils
