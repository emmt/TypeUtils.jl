module TypeUtils

export
    AbstractTypeStableFunction,
    TypeStableFunction,
    as,
    as_eltype,
    as_return,
    bare_type,
    convert_bare_type,
    convert_eltype,
    convert_floating_point_type,
    convert_real_type,
    destructure!,
    destructure,
    floating_point_type,
    parameterless,
    promote_eltype,
    real_type,
    restructure,
    return_type,
    struct_length,
    unitless

using Base: OneTo

if !isdefined(Base, :get_extension)
    using Requires
end

"""
    as(T, x)

yields `x` converted to type `T`.

"""
as(::Type{T}, x::T) where {T} = x
as(::Type{T}, x) where {T} = convert(T, x)::T

# Convert Cartesian index/indices to/from tuples.
as(::Type{Tuple}, x::CartesianIndex) = Tuple(x)
as(::Type{Tuple}, x::CartesianIndices) = x.indices
for X in (:CartesianIndex, :CartesianIndices)
    @eval begin
        # For more specific tuple types, first extract tuple contents, then
        # convert to the requested tuple type.
        as(::Type{T}, x::$X) where {T<:Tuple} = as(T, as(Tuple, x))

        # Use the constructors to convert tuples to Cartesian index/indices.
        as(::Type{$X}, x::Tuple) = $X(x)::$X
        as(::Type{$X{N}}, x::NTuple{N,Any}) where {N} = $X(x)::$X{N}
    end
end

# Conversion between symbols and strings is not supported by `convert`.
as(::Type{String}, x::Symbol) = String(x)
as(::Type{Symbol}, x::String) = Symbol(x)

"""
    f = as(T)

yields a callable object which converts its argument to type `T`. More
specifically, a call like `f(x)` yields `as(T, x)`.

"""
as(::Type{T}) where {T} = As{T}()

struct As{T} <: Function; end

(::As{T})(x) where {T} = as(T, x)

"""
    AbstractTypeStableFunction{T}

is the super-type of callable object with guaranteed returned type `T`.

"""
abstract type AbstractTypeStableFunction{T} <:Function end

struct TypeStableFunction{T,F} <: AbstractTypeStableFunction{T}
    callable::F
    TypeStableFunction{T}(f::F) where {T,F} = new{T,F}(f)
end

# Outer constructor.
function TypeStableFunction(f, argtypes::DataType...)
    T = Base.promote_op(f, argtypes...)
    return TypeStableFunction{T}(f)
end

# Conversion constructors.
TypeStableFunction{T}(f::TypeStableFunction{T}) where {T} = f
TypeStableFunction{T}(f::TypeStableFunction) where {T} = TypeStableFunction{T}(parent(f))

# Abstract constructor.
AbstractTypeStableFunction(f, argtypes::DataType...) = TypeStableFunction(f, argtypes...)
AbstractTypeStableFunction{T}(f) where {T} = TypeStableFunction{T}(f)

# Make instances of TypeStableFunction callable.
@inline (obj::TypeStableFunction{T})(args...; kwds...) where {T} =
    as(T, parent(obj)(args...; kwds...))

# Extend base methods.
Base.parent(obj::TypeStableFunction) = getfield(obj, :callable)

Base.return_types(::AbstractTypeStableFunction{T}; kwds...) where {T} = (T,)
Base.return_types(::AbstractTypeStableFunction{T}, ::DataType; kwds...) where {T} = (T,)

Base.promote_op(::AbstractTypeStableFunction{T}, ::DataType...) where {T} = T

for cls in (:AbstractTypeStableFunction, :TypeStableFunction,)
    @eval begin
        Base.convert(::Type{T}, f::T) where {T<:$(cls)} = f
        Base.convert(::Type{$(cls){T}}, f) where {T} = $(cls){T}(f)
    end
end

"""
    as_return(T, f) -> g
    TypeStableFunction{T}(f) -> g
    TypeStableFunction(f, argtypes...) -> g

yield a callable object `g` that wraps callable `f` for guaranteed returned
type `T`. Alternatively, the type(s) `argtypes...` of the function argument(s)
can be specified to infer the returned type `T`. Then the following holds:

    g(args...; kwds...) === as(T, f(args...; kwds...))

for any arguments `args...` and keywords `kwds...`. Note that due to limitation
of the `Base.promote_op` method, it is currently not possible to infer `T`
based on the types of the keywords.

 Methods [`return_type(g)`](@ref) and `parent(g)` may be used to retrieve `T`
and `f` respectively.

A similar object is given by:

    g = as(T)∘f

"""
as_return(::Type{T}, f) where {T} = TypeStableFunction{T}(f)

@doc @doc(as_return) TypeStableFunction

"""
    return_type(f, argtypes...) -> T

yields the type of the result returned by the callable object `f` when called
with arguments of types `argtypes...`.

See the warning in the documentation of `Base.promote_op` for the fragility of
such inference in some cases. There are no such issues if `f` is an instance of
[`TypeStableFunction`](@ref), e.g. built by [`as_return`][@ref), however
`argtypes...` are not checked for validity for such objects.

"""
return_type(::TypeStableFunction{T}, ::DataType...) where {T} = T
return_type(::Type{<:TypeStableFunction{T}}, ::DataType...) where {T} = T
return_type(f, argtypes::DataType...) = Base.promote_op(f, argtypes...)

"""
    parameterless(T)

yields the type `T` without parameter specifications. For example:

```julia
julia> parameterless(Vector{Float32})
Array
```

""" parameterless
#
# This subject is discussed in different places:
#
# - https://discourse.julialang.org/t/typeutils-dealing-with-types-in-julia/101584/13
#
# - https://stackoverflow.com/questions/42229901/getting-the-parameter-less-type
#
# The latter leads to define the method as (in old versions of Julia, the field
# name was `:primary`, but since Julia 0.7, it should be `:wrapper`):
#
#     @inline parameterless(::Type{T}) where {T} = getfield(Base.typename(T), :wrapper)
#
# The actual implementation is borrowed from `ConstructionBase` and should be
# less likely to be broken by internal changes in Julia:
#
parameterless(::Type{T}) where {T} = getfield(parentmodule(T), nameof(T))

"""
    TypeUtils.BareNumber

is the union of bare numeric types, that is `Real` or `Complex`.

"""
const BareNumber = Union{Real,Complex}

"""
    unitless(x)

yields `x` without its units if any. `x` may be a number or a numeric type. In
the latter case, `unitless` behaves like `bare_type`.

Compared to `ustrip` from the `Unitful` package, argument can be a numeric type
and, of course, `unitless` only requires the lightweight `TypeUtils` package to
be loaded.

"""
unitless(T::Type) = bare_type(T)
unitless(x::BareNumber) = x

"""
    bare_type(x) -> T <: Union{Real,Complex}

yields the bare numeric type `T` backing the storage of `x` which may be a
number or a numeric type. If `x` has units, they are discarded. Hence `T` is
always a unitless real or complex type.

Examples:

```jldoctest
julia> map(bare_type, (1, 3.14f0, π, 1 + 0im))
(Int64, Float32, Irrational{:π}, Complex{Int64})

julia> using Unitful

julia> map(bare_type, (u"3km/s", u"3.2km/s", typeof(u"2.1GHz")))
(Int64, Float64, Float64)
```

---
    bare_type(args...) -> T <: Union{Real,Complex}

yields the promoted bare numeric type of `args...`.

---
    bare_type() -> TypeUtils.BareNumber

yields the union of bare numeric types that may be returned by `bare_type` when
called with no arguments.

"""
bare_type() = BareNumber
bare_type(x::T) where {T} = bare_type(T)
bare_type(::Type{T}) where {T<:BareNumber} = T
bare_type(::Type{T}) where {T<:Number} = typeof(one(T))
@noinline bare_type(::Type{T}) where {T} =
    error("unknown bare numeric type for `", T, "`")

"""
    real_type(x) -> T <: Real

yields the bare numeric type `T` backing the storage of `x` which may be a
number of a numeric type. If `x` is a complex, `T` is the bare numeric type of
the real and imaginary parts of `x`. If `x` has units, they are discarded.
Hence `T` is always a unitless real type.

Examples:

```jldoctest
julia> using TypeUtils

julia> map(real_type, (-3.14f0, 1 + 0im, Complex{Int8}))
(Float32, Int64, Int8)

julia> using Unitful

julia> real_type(u"3km/s")
Int64
```

---
    real_type(args...)

yields the promoted bare real type of `args...`.

---
    real_type() -> Real

yields the supertype of the types that may be returned by `real_type` when
called with no arguments.

"""
real_type() = Real
real_type(x::T) where {T} = real_type(T)
real_type(::Type{T}) where {T<:Real} = T
real_type(::Type{Complex{T}}) where {T<:Real} = T
real_type(::Type{T}) where {T<:Number} = typeof(one(real(T)))
@noinline real_type(::Type{T}) where {T} = error("unknown bare real type for `", T, "`")

# Multiple arguments.
for f in (:bare_type, :real_type)
    @eval begin
        $f(a, b) = promote_type($f(a), $f(b))
        $f(a, b, c) = promote_type($f(a), $f(b), $f(c))
        @inline $f(a, b, c...) = promote_type($f(a), $f(b), map($f, c)...)
    end
end

"""
    convert_bare_type(T, x)

converts `x` so that its bare numeric type is that of `T`. Argument `x` may be
a number or a numeric type, while argument `T` must be a numeric type. If `x`
is one of `missing`, `nothing`, `undef`, or the type of one of these
singletons, `x` is returned.

This method may be extended with `T<:TypeUtils.BareNumber` and for `x` of
non-standard numeric type.

"""
convert_bare_type(::Type{T}, x) where {T<:Number} = convert_bare_type(bare_type(T), x)

# NOTE: All other specializations of `convert_bare_type(T,x)` are for `T<:BareNumber`.
convert_bare_type(::Type{T}, x::T) where {T<:BareNumber} = x
convert_bare_type(::Type{T}, x::BareNumber) where {T<:BareNumber} = convert(T, x)
@noinline convert_bare_type(::Type{T}, x) where {T<:BareNumber} = error(
   "unsupported conversion of bare numeric type of object of type `",
    typeof(x), "` to `", T, "`")

convert_bare_type(::Type{T}, ::Type{<:BareNumber}) where {T<:BareNumber} = T
@noinline convert_bare_type(::Type{T}, ::Type{S}) where {T<:BareNumber,S} = error(
    "unsupported conversion of bare numeric type of type `", S, "` to `", T, "`")

"""
    convert_real_type(T, x)

converts `x` so that its bare real type is that of `T`. Argument `x` may be a
number or a numeric type, while argument `T` must be a numeric type. If `x` is
one of `missing`, `nothing`, `undef`, or the type of one of these singletons,
`x` is returned.

This method may be extended with `T<:Real` and for `x` of non-standard numeric
type.

"""
convert_real_type(::Type{T}, x) where {T<:Number} = convert_real_type(real_type(T), x)

# NOTE: All other specializations of `convert_real_type(T,x)` are for `T<:Real`.
convert_real_type(::Type{T}, x::T) where {T<:Real} = x
convert_real_type(::Type{T}, x::Complex{T}) where {T<:Real} = x
convert_real_type(::Type{T}, x::Real) where {T<:Real} = convert(T, x)
convert_real_type(::Type{T}, x::Complex) where {T<:Real} = convert(Complex{T}, x)
@noinline convert_real_type(::Type{T}, x) where {T<:Real} = error(
    "unsupported conversion of bare real type of object of type `",
    typeof(x), "` to `", T, "`")

convert_real_type(::Type{T}, ::Type{<:Real}) where {T<:Real} = T
convert_real_type(::Type{T}, ::Type{<:Complex}) where {T<:Real} = Complex{T}
@noinline convert_real_type(::Type{T}, ::Type{S}) where {T<:Real,S} = error(
    "unsupported conversion of bare real type of type `", S, "` to `", T, "`")

# Special values/types.
const Special = Union{Missing,Nothing,typeof(undef)}
for (func, class) in ((:convert_bare_type, :BareNumber),
                      (:convert_real_type, :Real))
    @eval begin
        $func(::Type{<:$class}, x::Special) = x
        $func(::Type{<:$class}, ::Type{T}) where {T<:Special} = T
    end
end

"""
    floating_point_type(args...) -> T <: AbstractFloat

yields an appropriate floating-point type to represent the promoted numeric
type used by arguments `args...` for storing their value(s). Any units of the
arguments are ignored and the returned type is always unitless.

For numerical computations, a typical usage is:

    T = floating_point_type(x, y, ...)
    xp = convert_real_type(T, x)
    yp = convert_real_type(T, y)
    ...

to have numbers `x`, `y`, etc. converted to an appropriate common
floating-point type while preserving their units if any.

Also see [`real_type`](@ref) and [`convert_real_type`](@ref).

---
    floating_point_type() -> AbstractFloat

yields the supertype of the types that may be returned by `floating_point_type`
when called with no arguments.

"""
floating_point_type() = AbstractFloat
@inline floating_point_type(args...) = float(real_type(args...))

"""
    convert_floating_point_type(T, x)

converts `x` so that its bare real type is the floating-point type of `T`.
Argument `x` may be a number or a numeric type, while argument `T` must be a
numeric type. If `x` is one of `missing`, `nothing`, `undef`, or the type of
one of these singletons, `x` is returned.

This method may be extended with `T<:AbstractFloat` and for `x` of non-standard
numeric type.

"""
convert_floating_point_type(::Type{T}, x) where {T<:Number} =
    convert_real_type(floating_point_type(T), x)

"""
    promote_eltype(args...)

yields the promoted element type of its arguments. Arguments `args...` may be
anything implementing the `eltype` method.

"""
promote_eltype() = promote_type()
promote_eltype(arg) = eltype(arg)
@inline promote_eltype(args...) = promote_type(map(eltype, args)...)

"""
    convert_eltype(T, A) -> B

converts the elements of `A` to type `T`. The returned object is similar to `A`
except maybe for the element type. For example, if `A` is a range, then `B` is
also a range. If `T` is the element type of `A`, then `A` is returned.

Consider using [`as_eltype(T, A)`](@ref) to build an object that lazily
performs the conversion.

"""
convert_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractArray) where {T} = convert(AbstractArray{T}, A)

# Convert element type for tuples.
convert_eltype(::Type{T}, A::NTuple{N,T}) where {N,T} = A
convert_eltype(::Type{T}, A::Tuple) where {T} = map(as(T), A)

# Convert element type for Base.OneTo{T<:Integer} <: AbstractUnitRange{T}.
# Conversion to non-integer element types is handled by the more general method
# for AbstractUnitRange.
convert_eltype(::Type{T}, A::OneTo{T}) where {T<:Integer} = A
convert_eltype(::Type{T}, A::OneTo) where {T<:Integer} = OneTo{T}(last(A))

# Convert element type for AbstractUnitRange{T} <: OrdinalRange{T,T}.
convert_eltype(::Type{T}, A::AbstractUnitRange{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractUnitRange) where {T} =
    as(T, first(A)):as(T, last(A))

# Convert element type for other range types.
convert_eltype(::Type{T}, A::AbstractRange{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractRange) where {T} =
    as(T, first(A)):as(T, step(A)):as(T, last(A))

# Convert element type for LinRange{T,L<:Integer} <: AbstractRange{T}.
convert_eltype(::Type{T}, A::LinRange{T}) where {T} = A
convert_eltype(::Type{T}, A::LinRange) where {T} =
    LinRange(as(T, first(A)), as(T, last(A)), length(A))

"""
    as_eltype(T, A) -> B

yields an array which lazily converts its entries to type `T`. More
specifically, a call like `B[i]` yields `as(T,A[i])`.

Consider using [`convert_eltype(T, A)`](@ref) to perform the conversion once
and immediately.

"""
as_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
as_eltype(::Type{T}, A::AbstractArray) where {T} = AsEltype{T}(A)

struct AsEltype{T,N,L,A<:AbstractArray} <: AbstractArray{T,N}
    parent::A
end

AsEltype{T}(arr::A) where {T,N,A<:AbstractArray{<:Any,N}} =
    AsEltype{T,N,IndexStyle(A)===IndexLinear(),A}(arr)

Base.parent(A::AsEltype) = A.parent

# Implement abstract array API for `AsEltype` objects.
for func in (:axes, :length, :size)
    @eval Base.$func(A::AsEltype) = $func(parent(A))
end
Base.IndexStyle(::Type{<:AsEltype{<:Any,<:Any,true}}) = IndexLinear()
Base.IndexStyle(::Type{<:AsEltype{<:Any,<:Any,false}}) = IndexCartesian()

@inline function Base.getindex(A::AsEltype{T,N,true}, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds r = getindex(parent(A), i)
    return as(T, r)
end

@inline function Base.getindex(A::AsEltype{T,N,false}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds r = getindex(parent(A), I...)
    return as(T, r)
end

@inline function Base.setindex!(A::AsEltype{T,N,true}, x, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds setindex!(parent(A), x, i)
    return A
end

@inline function Base.setindex!(A::AsEltype{T,N,false}, x, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds setindex!(parent(A), x, I...)
    return A
end

"""
    destructure(obj) -> vals::Tuple

destructures object `obj` as a tuple of field values. Any structures in `obj`
are recursively destructured.

See also [`destructure!`](@ref), [`restructure`](@ref), and
[`struct_length`](@ref).

"""
function destructure end

@generated destructure(obj::T) where {T} = encode(destructure, :obj, T)

function encode(::typeof(destructure), base::Union{Symbol,Expr,QuoteNode}, ::Type{T}) where {T}
    expr = Expr(:tuple) # start with empty tuple
    encode!(destructure, expr.args, base, T)
    return quote
        $(Expr(:meta, :inline))
        return $expr
    end
end

function encode!(::typeof(destructure), code::AbstractVector,
                 base::Union{Symbol,Expr,QuoteNode}, ::Type{T}) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        for k in 1:fieldcount(T)
            encode!(destructure, code, :(getfield($base, $k)), fieldtype(T, k))
        end
    else
        push!(code, base)
    end
    nothing
end

"""
    destructure!(vals, obj; offset = firstindex(vals) - 1) -> vals

destructures object `obj` into `vals[offset+1:offset+n]` and returns `vals`. Here
`n = struct_length(obj)` is the total number of values stored by `obj`.

See also [`destructure`](@ref), [`restructure`](@ref), and
[`struct_length`](@ref).

"""
function destructure! end

@generated function destructure!(vals::AbstractVector, obj::T;
                                 offset::Integer = firstindex(vals) - 1) where {T}
    return encode(destructure!, :vals, :offset, :obj, T)
end

function encode(::typeof(destructure!), arr::Symbol, off::Symbol,
                base::Union{Symbol,Expr,QuoteNode}, ::Type{T}) where {T}
    code = Expr(:block, Expr(:meta, :inline))
    encode!(destructure!, code.args,
            (i, x) -> :($arr[$off + $i] = $x), # function to store each field
            base, T, 0)
    push!(code.args, :(return $arr))
    return code
end

function encode!(::typeof(destructure!), code::AbstractVector, f,
                 base::Union{Symbol,Expr,QuoteNode}, ::Type{T}, i::Int) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        for k in 1:fieldcount(T)
            i = encode!(destructure!, code, f, :(getfield($base, $k)), fieldtype(T, k), i)
        end
    else
        i += 1
        push!(code, f(i, base))
    end
    return i
end

"""
    restructure(T, vals; offset = firstindex(vals) - 1) -> obj::T

restructures values `vals[offset+1:offset+n]` into an object `obj` of type `T`.
Here `n = struct_length(T)` is the total number of values stored by an object
of type `T`.

The default constructors must exist for `T` and, recursively, for any
structured fields of `T`.

See also [`destructure`](@ref), [`destructure!`](@ref), and
[`struct_length`](@ref).

For an immutable concrete object `obj`, the following identity holds:

    restructure(typeof(obj), destructure(obj)) === obj

"""
function restructure end

@generated function restructure(::Type{T}, vals::Union{Tuple,AbstractVector};
                                offset::Integer = firstindex(vals) - 1) where {T}
    return encode(restructure, T, :vals, :offset)
end

function encode(::typeof(restructure), ::Type{T}, vals::Symbol, off::Symbol) where {T}
    code = Expr(:block, Expr(:meta, :inline))
    encode!(restructure, code, T, i -> :($vals[$off + $i]), 0)
    return code
end

function encode!(::typeof(restructure), code::Expr, ::Type{T}, f, i::Int) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        expr = if T <: Tuple
            Expr(:tuple)
        elseif isconcretetype(T)
            Expr(:call, T)
        else
            Expr(:call, parameterless(T))
        end
        for k in 1:fieldcount(T)
            i = encode!(restructure, expr, fieldtype(T, k), f, i)
        end
        if T <: Tuple
            push!(code.args, :(convert($T, $expr)))
        else
            push!(code.args, expr)
        end
    else
        i += 1
        push!(code.args, :($(f(i))))
    end
    return i
end

"""
    struct_length(T) -> n

yields the total number of values stored by the fields of a structured object
of type `T`. Argument may also be an object of type `T`.

"""
struct_length(x) = struct_length(typeof(x))

@generated function struct_length(::Type{T}) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        n = 0
        for k in 1:fieldcount(T)
            n += struct_length(fieldtype(T, k))
        end
        return n
    else
        return 1
    end
end

function __init__()
    @static if !isdefined(Base, :get_extension)
        # Extend methods to `Unitful` quantities when this package is loaded.
        @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" include(
            "../ext/TypeUtilsUnitfulExt.jl")
    end
end

end # module TypeUtils
