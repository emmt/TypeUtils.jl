module TypeUtils

export
    ArrayAxes,
    ArrayAxis,
    ArrayShape,
    RelaxedArrayShape,
    AbstractTypeStableFunction,
    TypeStableFunction,
    as,
    as_array_axes,
    as_array_axis,
    as_array_dim,
    as_array_shape,
    as_array_size,
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
    new_array,
    parameterless,
    promote_eltype,
    real_type,
    restructure,
    return_type,
    struct_length,
    to_same_type,
    to_same_concrete_type,
    unitless

using Base: OneTo

if !isdefined(Base, :get_extension)
    using Requires
end

"""
    @public args...

declares `args...` as being `public` even though they are not exported. For Julia version
< 1.11, this macro does nothing. Using this macro also avoid errors with CI and coverage
tools.

"""
macro public(args::Union{Symbol,Expr}...)
    VERSION ≥ v"1.11.0-DEV.469" ? esc(Expr(:public, args...)) : nothing
end

"""
    TypeUtils.Dim

is an alias to `eltype(Dims)`, the canonical integer type for an array dimension length.
In principle, `eltype(Dims) === Int` holds.

"""
const Dim = eltype(Dims)
@public Dim

"""
    ArrayAxis

is the canonical type of an array axis, an abstract unit range of `Int`s. Method
[`as_array_axis`](@ref) may be called to convert an argument to an array axis.

See also [`ArrayAxes`](@ref).

"""
const ArrayAxis = AbstractUnitRange{Dim}

"""
    ArrayAxes{N}

is the canonical type of array axes as returned by the base method
`axes(A::AbstractArray)`. It is an `N`-tuple of [`ArrayAxis`](@ref) instances. Method
[`as_array_axes`](@ref) may be called to convert arguments to array axes.

See also [`ArrayShape`](@ref), `Dims`.

"""
const ArrayAxes{N} = NTuple{N,ArrayAxis}

"""
    ArrayShape{N}

is the type of eligible argument to represent an `N`-dimensional array shape as usually
accepted in Julia: it is an `N`-tuple of integers and/or integer-valued unit ranges.

Methods [`as_array_shape`](@ref), [`as_array_axes`](@ref), and [`as_array_size`](@ref) may
be called on any `ArrayShape{N}` instance to convert it to a canonical form of array axes
or size.

Expression `eltype(ArrayShape)` yields the union of possible types for each entry of an
array shape. This may be used to specify a variable number of arguments possibly
representing an array shape.

See also [`ArrayAxes`](@ref), `Dims`.

"""
const ArrayShape{N} = NTuple{N,Union{Integer,AbstractUnitRange{<:Integer}}}

"""
    RelaxedArrayShape{N}

is like [`ArrayShape{N}`](@ref) but also includes `AbstractRange{<:Integer}}}`s (not just
`AbstractUnitRange{<:Integer}`s).

Methods [`as_array_shape`](@ref), [`as_array_axes`](@ref), and [`as_array_size`](@ref) may
be called on any argument of type `RelaxedArrayShape{N}` to convert it to a canonical form
of array axes or size. Likewise, methods [`as_array_dim`](@ref) and
[`as_array_axis`](@ref) may be called on any argument of type `eltype(RelaxedArrayShape)`
to convert it to a canonical array dimension length or axis. All these methods assert that
all ranges have a unit step.

!!! warning
    It is preferable to use [`ArrayShape{N}`](@ref) which better reflects Julia style.

"""
const RelaxedArrayShape{N} = NTuple{N,Union{Integer,AbstractRange{<:Integer}}}

"""
    TypeUtils.Unsupported(T::DataType...)

yields an union of types `T...` and of type `Unsupported`. Such an union can be used to
mark unsupported argument type(s) and yet define a method applicable with that type(s)
(presumably a method that throws an instructive error) and which can be extended later
with the same signature except that with `Unsupported(T...)` replaced by `T...`. This
trick avoids conflicts that prevent pre-compilation with package extensions.

For example, in the main package:

```julia
some_method(some_arg::Unsupported{SomeType}) =
    error("package `SomeOtherPackage` has not yet been loaded")
```

and in the extension (e.g. automatically loaded when using `SomeOtherPackage`):

```julia
some_method(some_arg::SomeType) = do_something_with(some_arg)
```

"""
struct Unsupported
    Unsupported() = error("it is not possible to instanciate this type")
end
Unsupported(T::DataType...) = Union{T...,Unsupported}
@public Unsupported

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
    R = T isa Union ? promote_type(ntuple(i -> getfield(T, i), Val(nfields(T)))...) : T
    isconcretetype(R) || throw(ArgumentError("cannot promote `$T` to a single concrete type"))
    return TypeStableFunction{R}(f)
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

yield a callable object `g` that wraps callable `f` for guaranteed returned type `T`.
Alternatively, the type(s) `argtypes...` of the function argument(s) can be specified to
infer the returned type `T`. Then the following holds:

    g(args...; kwds...) === as(T, f(args...; kwds...))

for any arguments `args...` and keywords `kwds...`. Note that due to limitation of the
`Base.promote_op` method, it is currently not possible to infer `T` based on the types of
the keywords.

Methods [`return_type(g)`](@ref) and `parent(g)` may be used to retrieve `T` and `f`
respectively.

If `T` is specified, a similar object is given by:

    g = as(T)∘f

"""
as_return(::Type{T}, f) where {T} = TypeStableFunction{T}(f)

@doc @doc(as_return) TypeStableFunction

"""
    as_array_shape(args...) -> r::Union{Dims,ArraysAxes}

converts array dimensions or ranges `args...` to a canonical form of array shape, one of:

* array size, that is a tuple of `Int`s. This is the result if all of `args...` are
  integers or instances of `Base.OneTo`, the latter, if any, being replaced by their
  lengths.

* array axes, that is a tuple of `AbstractUnitRange{Int}`s. This is the result if any of
  `args...` are non-`Base.OneTo` ranges, the integers being converted to instances of
  `Base.OneTo{eltype(Dims)}`.

The array dimensions or ranges may also be provided as a tuple.
[`RelaxedArrayShape{N}`](@ref) is the union of types of `N`-tuples to which
`as_array_shape` is applicable.

Also see [`as_array_size`](@ref), [`as_array_axes`](@ref), [`ArrayAxes`](@ref), `Dims`,
and [`new_array`](@ref).

"""
as_array_shape(::Tuple{}) = ()
as_array_shape(args::eltype(RelaxedArrayShape)...) = as_array_shape(args)
as_array_shape(args::Tuple{Vararg{Union{Integer,Base.OneTo{<:Integer}}}}) = as_array_size(args)
as_array_shape(args::RelaxedArrayShape) = as_array_axes(args)

"""
    as_array_size(args...) -> dims::Dims

converts array dimensions or ranges `args...` to a canonical form of array size, that is a
tuple of `eltype(Dims)`s. Any range in `args...` is replaced by its length.

The array dimensions or ranges may also be provided as a tuple.
[`RelaxedArrayShape{N}`](@ref) is the union of types of `N`-tuples to which
`as_array_size` is applicable.

Also see [`as_array_shape`](@ref), [`as_array_axes`](@ref), [`as_array_dim`](@ref),
`Dims`, and [`new_array`](@ref).

"""
as_array_size(::Tuple{}) = ()
as_array_size(args::eltype(RelaxedArrayShape)...) = as_array_size(args)
as_array_size(dims::Dims) = dims
as_array_size(args::RelaxedArrayShape) = map(as_array_dim, args)

"""
    as_array_axes(args...) -> rngs::ArrayAxes

converts array dimensions or ranges `args...` to a canonical form of array axes, that is a
tuple of `AbstractUnitRange{eltype(Dims)}`s. Any integer in `args...` is replaced by an
instance of `Base.OneTo{eltype(Dims)}`.

The array dimensions or ranges may also be provided as a tuple.
[`RelaxedArrayShape{N}`](@ref) is the union of types of `N`-tuples to which
`as_array_axes` is applicable.

Also see [`as_array_shape`](@ref), [`as_array_size`](@ref), [`as_array_axis`](@ref),
[`ArrayAxes`](@ref), `Dims`, and [`new_array`](@ref).

"""
as_array_axes(::Tuple{}) = ()
as_array_axes(args::eltype(RelaxedArrayShape)...) = as_array_axes(args)
as_array_axes(rngs::ArrayAxes) = rngs
as_array_axes(args::RelaxedArrayShape) = map(as_array_axis, args)

"""
    as_array_dim(arg) -> dim::eltype(Dims)

converts array dimension or range `arg` to a canonical array dimension, that is an
`eltype(Dims)`. If `arg` is a unit-step range, its length is returned.
[`eltype(RelaxedArrayShape)`](@ref RelaxedArrayShape) is the union of types to which
`as_array_dim` is applicable.

Also see [`as_array_size`](@ref), [`as_array_axis`](@ref), and `Dims`.

"""
as_array_dim(dim::Dim) = dim
as_array_dim(dim::Integer) = as(Dim, dim)
as_array_dim(rng::AbstractUnitRange{<:Integer}) = as_array_dim(length(rng))
as_array_dim(rng::AbstractRange{<:Integer}) =
     isone(step(rng)) ? as_array_dim(length(rng)) : throw_non_unit_step(rng)

"""
    as_array_axis(arg) -> rng::ArrayAxis

converts array dimension or range `arg` to a canonical array axis, that is an instance of
`AbstractUnitRange{eltype(Dims)}`. If `arg` is an integer, `Base.OneTo{eltype(Dims)}(arg)`
is returned. [`eltype(RelaxedArrayShape)`](@ref RelaxedArrayShape) is the union of types
to which `as_array_axis` is applicable.

Also see [`as_array_axes`](@ref), [`as_array_dim`](@ref), and `Dims`.

"""
as_array_axis(dim::Integer) = Base.OneTo{Dim}(dim)
as_array_axis(rng::ArrayAxis) = rng
as_array_axis(rng::AbstractUnitRange{<:Integer}) = ArrayAxis(rng)
as_array_axis(rng::AbstractRange{<:Integer}) =
     isone(step(rng)) ? UnitRange{Dim}(first(rng), last(rng)) : throw_non_unit_step(rng)

@noinline throw_non_unit_step(rng::AbstractRange) = throw(ArgumentError(
    "invalid non-unit step ($(step(rng))) for array axis"))

"""
    new_array(T, inds...) -> A
    new_array(T, (inds...,)) -> A

create a new array with element type `T` and shape defined by `inds...`, a list of array
dimension lengths and/or index ranges. The shape may also be specified as a tuple.

If `inds...` contains any index range other than `Base.OneTo`, an `OffsetArray{T}` is
returned; otherwise an `Array{T}` is returned. In the former case, an exception is thrown
if the package `OffsetArrays` has not been loaded.

Also see [`as_array_shape`](@ref), [`as_array_axes`](@ref), and [`as_array_size`](@ref).

"""
new_array(::Type{T}, inds::eltype(RelaxedArrayShape)...) where {T} = new_array(T, inds)
new_array(::Type{T}, inds::RelaxedArrayShape) where {T} = new_array(T, as_array_shape(inds))
new_array(::Type{T}, dims::Dims{N}) where {T,N} = Array{T,N}(undef, dims)
new_array(::Type{T}, rngs::Unsupported(ArrayAxes{N})) where {T,N} =
    error("package `OffsetArrays` must be loaded for such array index ranges")

# The following is needed to yield regular arrays if possible.
new_array(::Type{T}, inds::NTuple{N,Union{Integer,Base.OneTo}}) where {T,N} =
    new_array(T, as_array_size(inds))

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
    to_same_type(x1, x2, ...) -> xp1, xp2, ...

converts instances `x1`, `x2`, ... to the same type. This method may be used instead of
`promote(x1,x2,...)` which does not warrant that the converted values have the same type.

Example:

```julia
julia> using Unitful, TypeUtils

julia> promote(2, 3.0)
(2.0, 3.0)

julia> to_same_type(2, 3.0)
(2.0, 3.0)

julia> promote(2u"mm", 3.0)
(2.0 mm, 3.0)

julia> to_same_type(2u"mm", 3.0)
ERROR: ArgumentError: types `Quantity{Int64, 𝐋, Unitful.FreeUnits{(mm,), 𝐋, nothing}}` and `Float64` cannot be converted to a common concrete type
Stacktrace:
 ...

julia> promote(2u"mm", 4.0u"nm")
(0.002 m, 4.0e-9 m)

julia> to_same_type(2u"mm", 4.0u"nm")
(0.002 m, 4.0e-9 m)
```

Also see [`to_same_concrete_type`](@ref).

""" function to_same_type end

# No conversion needed for arguments of the same type.
to_same_type() = ()
to_same_type(x) = (x,)
to_same_type(x1::T, x2::T) where {T} = (x1, x2)
to_same_type(xs::T...) where {T} = xs

# For arguments of different types, the promoted type must be concrete otherwise it cannot
# have direct instances.
function to_same_type(x1::T1, x2::T2) where {T1,T2}
    T = to_same_concrete_type(T1, T2)
    return as(T, x1), as(T, x2)
end
@inline function to_same_type(xs...)
    T = to_same_concrete_type(map(typeof, xs)...)
    return map(as(T), xs)
end

## Error catcher, arguments must be instances not types.
to_same_type(xs::DataType...) =
    throw(ArgumentError("argument(s) must be instance(s) not type(s)"))

"""
    to_same_concrete_type(T1::Type, T2::Type, ...) -> T::Type

yields `T = promote_type(T1, T2, ...)` throwing an exception if `T` is not a concrete
type.

Also see [`to_same_type`](@ref).

""" function to_same_concrete_type end

to_same_concrete_type() = throw(ArgumentError("no type(s) specified"))

function to_same_concrete_type(::Type{T}) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

function to_same_concrete_type(::Type{T}, ::Type{T}) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

function to_same_concrete_type(::Type{T1}, ::Type{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    isconcretetype(T) || throw_no_common_concrete_type(T1, T2)
    return T
end

@inline function to_same_concrete_type(::Type{T}...) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

@inline function to_same_concrete_type(Ts::Type...)
    T = promote_type(Ts...)
    isconcretetype(T) || throw_no_common_concrete_type(Ts...)
    return T
end

@noinline throw_not_concrete_type(T::Type) =
    throw(ArgumentError("type `$T` is not a concrete type"))

@noinline throw_no_common_concrete_type(T1::Type, T2::Type) =
    throw(ArgumentError("types `$T1` and `$T2` cannot be converted to a common concrete type"))

@noinline throw_no_common_concrete_type(Ts::Type...) =
    throw(ArgumentError(*("types `", join(Ts, "`, `", "`, and `"),
                          "` cannot be converted to a common concrete type")))

"""
    TypeUtils.BareNumber

is the union of bare numeric types, that is `Real` or `Complex`.

"""
const BareNumber = Union{Real,Complex}
@public BareNumber

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

converts the element type of object/type `A` to type `T`. The returned object/type is
similar to `A` except maybe for the element type. For example, if `A` is a range, then `B`
is also a range. If `T` is the element type of `A`, then `A` may be returned.

Consider using [`as_eltype(T, A)`](@ref) to build an object that lazily performs the
conversion.

To simplify extending `convert_eltype` for objects `A` of given type, the default behavior
is:

    convert_eltype(T, A) = as(convert_eltype(T, typeof(A)), A)

so that it may be sufficient to extend `convert_eltype` for the type of the objects.

"""
convert_eltype(::Type{T}, x::X) where {T,X} = as(convert_eltype(T, X), x)
convert_eltype(::Type{T}, ::Type{X}) where {T,X} =
    error("don't know how to convert the element type of type `$X` to `$T`")

convert_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractArray) where {T} = as(AbstractArray{T}, A)

# Convert element type for numbers.
convert_eltype(::Type{T}, ::Type{<:Number}) where {T} = T

# Convert element type for tuples. See `_countuple` in `base/tuple.jl` for the best
# way to extract the number of elements in a tuple given its type.
convert_eltype(::Type{T}, ::Type{<:NTuple{N,Any}}) where {N,T} = NTuple{N,T}
convert_eltype(::Type{T}, A::NTuple{N,T}) where {N,T} = A
convert_eltype(::Type{T}, A::Tuple) where {T} = map(as(T), A)

# Convert element type for `Base.OneTo` and `UnitRange`. For `T` non-integer,
# a `Base.OneTo` instance becomes a `UnitRange` one (in a predictible way).
convert_eltype(::Type{T}, ::Type{<:OneTo}) where {T<:Integer} = OneTo{T}
convert_eltype(::Type{T}, ::Type{<:OneTo}) where {T} = UnitRange{T}
convert_eltype(::Type{T}, ::Type{<:UnitRange}) where {T} = UnitRange{T}
convert_eltype(::Type{T}, A::Union{OneTo,UnitRange}) where {T} =
    as(convert_eltype(T, typeof(A)), A)

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
        # Extend methods when other packages are loaded.
        @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" include(
            "../ext/TypeUtilsUnitfulExt.jl")
        @require OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881" include(
            "../ext/TypeUtilsOffsetArraysExt.jl")
    end
end

end # module TypeUtils
