"""
    get_precision(x) -> T<:AbstractFloat
    get_precision(typeof(x)) -> T<:AbstractFloat

yield the numerical precision of number/object `x`. If `x` is a floating-point value, its
floating-point type is returned; if `x` stores floating-point values, their promoted
floating-point type is returned; otherwise, `AbstractFloat` is returned.

For type-stability, the precision shall be a *trait* which does only depend of the type of
`x`.

Basically, `get_precision` supports numbers, and arrays or tuples of numbers. It can be
specialized for other object types defined in foreign packages. For objects of type, say
`SomeType`, it is only needed to extend the method:

```julia
TypeUtils.get_precision(::Type{T}) where {T<:SomeType} = ...
```

See also [`adapt_precision`](@ref) and [`TypeUtils.Precision`](@ref).

"""
get_precision(x::T) where {T} = get_precision(T)
get_precision(::Type) = AbstractFloat # pass-through
get_precision(::Type{T}) where {T<:Precision} = T
get_precision(::Type{<:Complex{T}}) where {T} = get_precision(T)
get_precision(::Type{<:AbstractArray{T}}) where {T} = get_precision(T)
get_precision(::Type{<:Factorization{T}}) where {T} = get_precision(T)

@generated function get_precision(::Type{T}) where {T<:Union{Tuple,NamedTuple}}
    # NOTE Using a `Ref` for `r` or `t` here is significantly slower.
    r = AbstractFloat
    for s in T.types
        t = get_precision(s)::Type{<:AbstractFloat}
        if isconcretetype(t)
            if r == AbstractFloat
                r = t
            else
                r = promote_type(r, t)
            end
        end
    end
    return r
end

"""
    adapt_precision(T::Type{<:AbstractFloat}, x) -> y

yields an object `y` similar to `x` but with numerical precision specified by the
floating-point type `T`. If `x` has already the required precision or if setting its
precision is irrelevant or not implemented, `x` is returned unchanged. Setting the
precision shall not change the dimensions of dimensionful numbers. If `T` is
`AbstractFloat`, the default floating-point type [`TypeUtils.default_precision`](@ref) is
assumed.

For a number `x`, `adapt_precision(T, x)` behaves as [`convert_real_type(T, x)`](@ref
convert_real_type) and `adapt_precision(T, typeof(x))` may be used to infer the
corresponding type with precision `T`.

Example:

```jldoctest; setup=:(using TypeUtils)
julia> adapt_precision(Float32, (1, 0x07, ("hello", 1.0, 3.0 - 2.0im, π)))
(1.0f0, 7.0f0, ("hello", 1.0f0, 3.0f0 - 2.0f0im, 3.1415927f0))
```

As can be seen, all numerical values are converted.

Basically, `adapt_precision` supports numbers, and arrays or tuples of numbers. It can be
specialized for other object types defined in foreign packages by specializing:

```julia
TypeUtils.adapt_precision(::Type{T}, x::SomeType) where {T<:TypeUtils.Precision} = ...
```

where `SomeType` is the object type and where the restriction `T<:TypeUtils.Precision` is
to make sure the above method is only called with a concrete floating-point type `T`.

See also [`get_precision`](@ref), [`convert_real_type`](@ref),
[`TypeUtils.Precision`](@ref), and [`TypeUtils.default_precision`](@ref).

"""
adapt_precision(::Type{AbstractFloat}, x::Any) = adapt_precision(default_precision, x)
adapt_precision(::Type{T}, x::Any) where {T} = throw_not_precision(T)

@noinline throw_not_precision(::Type{T}) where {T} = throw(ArgumentError(
    "precision `$T` is not a concrete floating-point type"))

# Pass-through by default.
adapt_precision(::Type{T}, x::Any) where{T<:Precision} = x
adapt_precision(::Type{T}, ::Type{S}) where {T<:Precision,S<:Any} = S

# For bare numbers and bare numerical types, `adapt_precision` behaves like
# `convert_real_type`.
adapt_precision(::Type{T}, x::T) where {T<:Precision} = x
adapt_precision(::Type{T}, x::BareNumber) where {T<:Precision} = convert_real_type(T, x)
adapt_precision(::Type{T}, ::Type{T}) where {T<:Precision} = T
adapt_precision(::Type{T}, ::Type{S}) where {T<:Precision,S<:BareNumber} =
    convert_real_type(T, S)

# Adapt precision of array types.
adapt_precision(::Type{T}, ::Type{A}) where {T<:Precision,A<:AbstractArray{T}} = A
adapt_precision(::Type{T}, ::Type{A}) where {T<:Precision,A<:AbstractArray} =
    convert_eltype(adapt_precision(T, eltype(A)), A)

# Adapt precision of numerical arrays.
adapt_precision(::Type{T}, A::AbstractArray{T}) where {T<:Precision} = A
adapt_precision(::Type{T}, A::AbstractArray{S}) where {T<:Precision,S} =
    convert_eltype(adapt_precision(T, S), A)

# Adapt precision of factorizations.
adapt_precision(::Type{T}, A::Factorization{T}) where {T<:Precision} = A
adapt_precision(::Type{T}, A::Factorization{S}) where {T<:Precision,S} =
    convert_eltype(adapt_precision(T, S), A)

# Set precision for tuples.
adapt_precision(::Type{T}, x::Union{Tuple,NamedTuple}) where {T<:Precision} =
    map(adapt_precision(T), x)

"""
    f = adapt_precision(T)

builds a callable object `f` such that `f(x)` is equivalent to `adapt_precision(T, x)`. If
`T` is `AbstractFloat`, the default floating-point type
[`TypeUtils.default_precision`](@ref) is assumed.

"""
adapt_precision(::Type{AbstractFloat}) = adapt_precision(default_precision)
adapt_precision(::Type{T}) where {T<:Precision} = Converter(adapt_precision, T)
adapt_precision(::Type{T}) where {T} = throw_not_precision(T)
