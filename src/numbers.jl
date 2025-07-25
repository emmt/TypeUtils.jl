# Methods related to numbers.

"""
    bare_type(x) -> T <: Union{Real,Complex}

yields the bare numeric type `T` backing the storage of `x` which may be a number or a
numeric type. If `x` has units, they are discarded. Hence `T` is always a dimensionless real
or complex type.

Examples:

```jldoctest; setup=:(using TypeUtils)
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

yields the union of bare numeric types that may be returned by `bare_type` when called with
no arguments.

"""
bare_type() = BareNumber
bare_type(x::T) where {T} = bare_type(T)
bare_type(::Type{T}) where {T<:BareNumber} = T
bare_type(::Type{T}) where {T<:Number} = typeof(one(T))
@noinline bare_type(::Type{T}) where {T} =
    error("unknown bare numeric type for `", T, "`")

"""
    real_type(x) -> T <: Real

yields the bare numeric type `T` backing the storage of `x` which may be a number of a
numeric type. If `x` is a complex, `T` is the bare numeric type of the real and imaginary
parts of `x`. If `x` has units, they are discarded. Hence `T` is always a dimensionless real
type.

Examples:

```jldoctest; setup=:(using TypeUtils)
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

yields the supertype of the types that may be returned by `real_type` when called with no
arguments.

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
        @inline $f(a, b, c...) = promote_type($f(a, b), $f(c...))
    end
end

"""
    convert_bare_type(T, x)

converts `x` so that its bare numeric type is that of `T`. Argument `x` may be a number or a
numeric type, while argument `T` must be a numeric type. If `x` is one of `missing`,
`nothing`, `undef`, or the type of one of these singletons, `x` is returned.

This method may be extended with `T<:TypeUtils.BareNumber` and for `x` of non-standard
numeric type.

"""
convert_bare_type(::Type{T}, x) where {T} = convert_bare_type(bare_type(T), x)

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
    convert_bare_type(T) -> f

yields a callable object `f` such that `f(x)` yields `convert_bare_type(T, x)` for any
`x`.

"""
convert_bare_type(::Type{T}) where {T} = Converter(convert_bare_type, bare_type(T))

"""
    convert_real_type(T, x)

converts `x` so that its bare real type is that of `T`. Argument `x` may be a number or a
numeric type, while argument `T` must be a numeric type. If `x` is one of `missing`,
`nothing`, `undef`, or the type of one of these singletons, `x` is returned.

This method may be extended with `T<:Real` and for `x` of non-standard numeric type.

"""
convert_real_type(::Type{T}, x) where {T} = convert_real_type(real_type(T), x)

# NOTE: All other specializations of `convert_real_type(T,x)` are for `T<:Real`.
convert_real_type(::Type{T}, x::T) where {T<:Real} = x
convert_real_type(::Type{T}, x::Complex{T}) where {T<:Real} = x
convert_real_type(::Type{T}, x::Real) where {T<:Real} = as(T, x)
convert_real_type(::Type{T}, x::Complex) where {T<:Real} = as(Complex{T}, x)
@noinline convert_real_type(::Type{T}, x) where {T<:Real} = error(
    "unsupported conversion of bare real type of object of type `",
    typeof(x), "` to `", T, "`")

convert_real_type(::Type{T}, ::Type{<:Real}) where {T<:Real} = T
convert_real_type(::Type{T}, ::Type{<:Complex}) where {T<:Real} = Complex{T}
@noinline convert_real_type(::Type{T}, ::Type{S}) where {T<:Real,S} = error(
    "unsupported conversion of bare real type of type `", S, "` to `", T, "`")

"""
    convert_real_type(T) -> f

yields a callable object `f` such that `f(x)` yields `convert_real_type(T, x)` for any
`x`.

"""
convert_real_type(::Type{T}) where {T} = Converter(convert_real_type, real_type(T))

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

yields an appropriate floating-point type to represent the promoted numeric type of
arguments `args...` for storing their value(s). Any units of the arguments are ignored and
the returned type is always dimensionless.

For numerical computations, a typical usage is:

    T = floating_point_type(x, y, ...)
    xp = convert_real_type(T, x)
    yp = convert_real_type(T, y)
    ...

to have numbers `x`, `y`, etc. converted to an appropriate common floating-point type
while preserving their units if any.

Also see [`real_type`](@ref) and [`convert_real_type`](@ref).

---
    floating_point_type() -> AbstractFloat

yields the supertype of the types that may be returned by `floating_point_type` when
called with no arguments.

"""
floating_point_type() = AbstractFloat
@inline floating_point_type(args...) = float(real_type(args...))

"""
    convert_floating_point_type(T, x)

converts `x` so that its bare real type is the floating-point type of `T`. Argument `x` may
be a number or a numeric type, while argument `T` must be a numeric type. If `x` is one of
`missing`, `nothing`, `undef`, or the type of one of these singletons, `x` is returned.

This method may be extended with `T<:AbstractFloat` and for `x` of non-standard numeric
type.

"""
convert_floating_point_type(::Type{T}, x) where {T} =
    convert_real_type(floating_point_type(T), x)

"""
    convert_floating_point_type(T) -> f

yields a callable object `f` such that `f(x)` yields `convert_floating_point_type(T, x)` for
any `x`.

"""
convert_floating_point_type(::Type{T}) where {T} =
    Converter(convert_floating_point_type, floating_point_type(T))

"""
    assert_floating_point(Bool, x) -> bool

yields whether `x` uses floating-point to store its value(s). For n-tuples, the same
floating-point type must be used for all values.

"""
assert_floating_point(::Type{Bool}, x) =
    assert_floating_point(Bool, typeof(x))

assert_floating_point(::Type{Bool}, ::Type{T}) where {T<:Union{Tuple,AbstractArray}} =
    assert_floating_point(Bool, eltype(T))

assert_floating_point(::Type{Bool}, ::DataType) = false

assert_floating_point(::Type{Bool}, ::Type{T}) where {T<:Number} =
    _assert_floating_point(Bool, real_type(T))

_assert_floating_point(::Type{Bool}, ::Type{T}) where {T<:AbstractFloat} = isconcretetype(T)
_assert_floating_point(::Type{Bool}, ::Type{T}) where {T<:Number} = false

"""
    assert_floating_point([name,] x)

throws an exception if `x` does not use floating-point to store its value(s). Optional
`name` argument is to specify the name to use for the error message.

See also [`@assert_floating_point`](@ref).

"""
assert_floating_point(x) = assert_floating_point(typeof(x))
assert_floating_point(::Type{T}) where {T} =
    assert_floating_point(Bool, T) ? nothing : throw_not_floating_point(T)

assert_floating_point(name::Union{Symbol,AbstractString}, x) = assert_floating_point(name, typeof(x))
assert_floating_point(name::Union{Symbol,AbstractString}, ::Type{T}) where {T} =
    assert_floating_point(Bool, T) ? nothing : throw_not_floating_point(name, T)

@noinline throw_not_floating_point(::Type{T}) where {T} =
    throw(ArgumentError("type `$T` is not floating-point"))
@noinline throw_not_floating_point(name::Union{Symbol,AbstractString}, ::Type{T}) where {T} =
    throw(ArgumentError("argument or variable `$name` of type `$T` is not floating-point"))

"""
    nearest(T::Type, x) -> y::T

yields the value or instance of type `T` that is the nearest to `x`. For `T` integer and `x`
real, it can be seen as rounding with clamping.

"""
nearest(::Type{T}, x::T) where {T} = x
nearest(::Type{T}, x::Real) where {T<:Real} = as(T, x) # by default, simply convert...

# Generic methods for integer to nearest integer.
nearest(::Type{T}, x::T) where {T<:Integer} = x
function nearest(::Type{T}, x::Integer) where {T<:Integer}
    # Use `ifelse` here because conversion by `%` does not throw.
    lo, hi = typemin(T), typemax(T)
    return ifelse(x ≤ lo, lo, ifelse(x ≥ hi, hi, x % T))
end

# Optimized methods for integer to nearest integer.
for T in BIT_INTEGERS, X in BIT_INTEGERS
    if X === T
        # No conversion needed.
        @eval nearest(::Type{$T}, x::$T) = x
    elseif typemin(T) ≤ typemin(X) && typemax(X) ≤ typemax(T)
        # No clamping need, just convert with `rem`.
        @eval nearest(::Type{$T}, x::$X) = x % $T
    end
end

# Non-integer real to nearest integer.
function nearest(::Type{T}, x::Real) where {T<:Integer}
    # We cannot use `ifelse` here because conversion may throw an `InexactError`. This will
    # occur anyway for `NaN`s.
    r = round(x)
    if r ≤ (lo = typemin(T))
        return lo
    elseif r ≥ (hi = typemax(T))
        return hi
    else
        return T(r)::T
    end
end

# Special case of `Bool`.
nearest(::Type{Bool}, x::Integer) = x > zero(x)
function nearest(::Type{Bool}, x::Real)
    isnan(x) && throw(InexactError(:nearest, Bool, x))
    r = round(x)
    return r > zero(r)
end

# Special case of `BigInt` for which `typemin` and `typemax` make no sense.
nearest(::Type{BigInt}, x::BigInt) = x
nearest(::Type{BigInt}, x::Integer) = BigInt(x)
nearest(::Type{BigInt}, x::Real) = round(BigInt, x)
nearest(::Type{BigInt}, x::Irrational) = round(BigInt, round(x))

"""
    nearest(T::Type) -> f

yields a callable object `f` such that `f(x)` yields `nearest(T, x)` for any `x`.

"""
nearest(::Type{T}) where {T} = Converter(nearest, T)

"""
    unitless(x)

yields `x` without its units if any. `x` may be a number or a numeric type. In the latter
case, `unitless` behaves like `bare_type`.

Compared to `ustrip` from the `Unitful` package, argument can be a numeric type and, of
course, `unitless` only requires the lightweight `TypeUtils` package to be loaded.

"""
unitless(T::Type) = bare_type(T)
unitless(x::BareNumber) = x
