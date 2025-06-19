"""
    TypeUtils.Dim

is an alias to `eltype(Dims)`, the canonical integer type for an array dimension length.
In principle, `eltype(Dims) === Int` holds.

"""
const Dim = eltype(Dims)

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

# Type of array shape that yields a regular `Array` in `new_array`.
const RegularArrayShape{N} = NTuple{N,Union{Integer,Base.OneTo{<:Integer}}}

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
    using TypeUtils: Unsupported

    Unsupported(T)
    Union{Unsupported, T}
    Unsupported(T1, T2, ...)
    Union{Unsupported, T1, T2, ...}

yield an union that can be used to define a method applicable with unsupported argument
type(s) `T` or `T1`, `T2`, and `...` (presumably a method that throws an instructive
error) and which can be extended later with the same signature except that with
`Union{Unsupported, T}` or `Unsupported(T)` replaced by `T` and
`Union{Unsupported,T1,T2,...}` or `Unsupported(T1,T2,...)` replaced by `Union{T1,T2,...}`.
This trick avoids conflicts that prevent pre-compilation with package extensions.

If `T` or any of `T1`, `T2`, and `...` have unspecified type parameters, the
`Union{Unsupported,...}` syntax may have to be used.

For example, in the main package:

```julia
using TypeUtils: Unsupported
some_method(some_arg::Union{Unsupported,SomeType}) =
    error("package `SomeOtherPackage` has not yet been loaded")
```

and in the extension (e.g. automatically loaded when using `SomeOtherPackage`):

```julia
some_method(some_arg::SomeType) = do_something_with(some_arg)
```

"""
struct Unsupported
    Unsupported() = error("it is not possible to instantiate this type")
end
Unsupported(T::DataType...) = Union{Unsupported,T...}

# Lists of machine integer types.
function _bit_integers(pfx)
    # FIXME This is equivalent to:
    #
    # map(eval, Iterators.takewhile(Base.Fix1(isdefined, Base),
    #                               Iterators.map(i -> Symbol(pfx, 8<<i),
    #                                             Iterators.countfrom(0))))
    #
    # but not all these iterators exist in old Julia versions.
    types = Type[]
    nbits = 8
    while true
        sym = Symbol(pfx, nbits)
        isdefined(Base, sym) || return types
        push!(types, eval(sym))
        nbits *= 2
    end
end
const BIT_INTEGERS = (Bool, _bit_integers("Int")..., _bit_integers("UInt")...,)
# FIXME `filter` does not work on collections other than arrays and dictionaries for
#       Julia < 1.4. Hence, `collect`...
const SIGNED_BIT_INTEGERS = (filter(T -> T <: Signed, collect(BIT_INTEGERS))...,)
const UNSIGNED_BIT_INTEGERS = (filter(T -> T <: Unsigned, collect(BIT_INTEGERS))...,)

"""
    TypeUtils.BareNumber

is the union of bare numeric types, that is `Real` or `Complex`.

"""
const BareNumber = Union{Real,Complex}

"""
    Precision

is the union of concrete floating-point types that can be used to specify the precision
with [`adapt_precision`](@ref).

"""
const Precision = Union{Float16,Float32,Float64,BigFloat}

const default_precision = Float64

"""
    c = TypeUtils.Converter(f, T::Type)

builds a lightweight callable object `c` such that `c(x)` yields `f(T, x)` for any `x`.
Converter objects are suitable to map conversion to type `T` by function `f` to
collections.

This is similar to `Base.Fix1(f,T)` except that `sizeof(Base.Fix1(f,T)) = sizeof(Int)`
while `sizeof(TypeUtils.Converter(f,T)) = 0`.

"""
struct Converter{F,T} <: Function
    func::F
    Converter(func::F, ::Type{T}) where {F,T} = new{F,T}(func)
end
(c::Converter{F,T})(x) where {F,T} = c.func(T, x)

struct AsEltype{T,N,L,A<:AbstractArray} <: AbstractArray{T,N}
    parent::A
    AsEltype{T}(arr::A) where {T,N,A<:AbstractArray{<:Any,N}} =
        new{T,N,IndexStyle(A)===IndexLinear(),A}(arr)
end

"""
    AbstractTypeStableFunction{T}

is the super-type of callable object with guaranteed returned type `T`.

"""
abstract type AbstractTypeStableFunction{T} <:Function end

struct TypeStableFunction{T,F} <: AbstractTypeStableFunction{T}
    callable::F
    TypeStableFunction{T}(f::F) where {T,F} = new{T,F}(f)
end
