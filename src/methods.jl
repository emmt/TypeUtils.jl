# General conversion methods.

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
        # For more specific tuple types, first extract tuple contents, then convert to the
        # requested tuple type.
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
    as(T) -> f

yields a callable object `f` such that `f(x)` yields `as(T, x)` for any `x`.

"""
as(::Type{T}) where {T} = Converter(as, T)

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
    to_same_concrete_type(Ts::Type...) -> T::Type

yields `T = promote_type(Ts...)` throwing an exception if `T` is not a concrete type.

Also see [`to_same_type`](@ref).

""" function to_same_concrete_type end

to_same_concrete_type() = throw(ArgumentError("no type(s) specified"))

@inline function to_same_concrete_type(::Type{T}...) where {T}
    isconcretetype(T) || throw_not_concrete_type(T)
    return T
end

function to_same_concrete_type(::Type{T1}, ::Type{T2}) where {T1,T2}
    T = promote_type(T1, T2)
    isconcretetype(T) || throw_no_common_concrete_type(T1, T2)
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
# The latter leads to define the method as (in old versions of Julia, the field name was
# `:primary`, but since Julia 0.7, it should be `:wrapper`):
#
#     @inline parameterless(::Type{T}) where {T} = getfield(Base.typename(T), :wrapper)
#
# The actual implementation is borrowed from `ConstructionBase` and should be less likely to
# be broken by internal changes in Julia:
#
parameterless(::Type{T}) where {T} = getfield(parentmodule(T), nameof(T))

"""
    is_signed(x)
    is_signed(typeof(x))

yield whether number `x` can be negated while retaining the same type. The result only
depends on the type of `x`. The result is `false` for `typeof(x) <:
Union{U,Rational{U},Complex{U}} where {U<:Union{Bool,Unsigned}}` as well as quantities based
on these bare numeric types.

"""
is_signed(x::Number) = is_signed(typeof(x))
is_signed(::Type{T}) where {T} = _is_signed(bare_type(T))
_is_signed(::Type{<:Number}) = true
_is_signed(::Type{<:Union{Bool,Unsigned}}) = false
_is_signed(::Type{<:Rational{T}}) where {T} = is_signed(T)
_is_signed(::Type{<:Complex{T}}) where {T} = is_signed(T)
