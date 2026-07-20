"""
    get_precision(x) -> T<:AbstractFloat
    get_precision(typeof(x)) -> T<:AbstractFloat

Return the numerical precision of number/object `x`. If `x` is a floating-point value, its
floating-point type is returned; if `x` stores floating-point values, their promoted
floating-point type is returned; otherwise, `AbstractFloat` is returned.

By default, `get_precision` is applicable to any object or type and assuming that the
precision is a *trait*, hence, which only depends on the type of `x`. This ensures
*type-stability*. However, if a specific behavior must be implemented to a type, say,
`SomeType`, it is only needed to extend the method:

```julia
TypeUtils.get_precision(::Type{T}) where {T<:SomeType} = ...
```

which shall return one of the floating-point types of union `TypeUtils.Precision` or
`AbstractFloat`.

If 0, 2 or more arguments are specified, `get_precision` returns the promoted precision of
all given arguments.

See also [`adapt_precision`](@ref) and [`TypeUtils.Precision`](@ref).

"""
get_precision(x::Any) = get_precision(typeof(x)) # precision is a trait
get_precision(::Type{T}) where {T<:Precision} = T
get_precision(::Type{T}) where {T<:Real} = AbstractFloat
get_precision(::Type{T}) where {T<:TypeVar} = AbstractFloat

# The following specialization is needed for abstract arrays like `StepRangeLen` which have
# fields of higher precision to avoid rounding errors.
get_precision(::Type{<:AbstractArray{T}}) where {T} = get_precision(T)

# Get precision for 0 or more than 2 arguments.
get_precision() = AbstractFloat
@inline get_precision(x, y, z...) = get_precision(get_precision(x, y), z...)

# Get precision of exactly 2 types. Type assertion in the most generic method is to avoid
# stack-overflow for an invalid specialization.
get_precision(::Type{AbstractFloat}, ::Type{AbstractFloat}) = AbstractFloat
get_precision(::Type{T}, ::Type{AbstractFloat}) where {T<:AbstractFloat} = get_precision(T)
get_precision(::Type{AbstractFloat}, ::Type{T}) where {T<:AbstractFloat} = get_precision(T)
get_precision(::Type{S}, ::Type{T}) where {S<:AbstractFloat,T<:AbstractFloat} =
    get_precision(promote_type(S, T))
get_precision(x, y) = get_precision(get_precision(x)::Type{<:AbstractFloat},
                                    get_precision(y)::Type{<:AbstractFloat})

# Return a simple vector of type parameters. (Must be defined before any generated function
# that may use it.)
function type_parameters(::Type{T}) where {T}
    try
        getfield(T, :parameters)
    catch
        try
            getfield(getfield(T, :body), :parameters)
        catch
            Base.Core.svec()
        end
    end
end

# Default method for structure, tuple, and named tuple types.
@generated function get_precision(::Type{S}) where {S}
    T = AbstractFloat
    for p in type_parameters(S)
        T = get_precision(T, p)
    end
    return T
end

"""
    adapt_precision(T::Type{<:AbstractFloat}, x) -> y

Return an object `y` similar to `x` but with numerical precision specified by the
floating-point type `T`. If `x` has already the required precision or if setting its
precision is irrelevant or not implemented, `x` is returned unchanged. Setting the precision
shall not change the units of dimensionful numbers. If `T` is `AbstractFloat`, the default
floating-point type [`TypeUtils.default_precision`](@ref) is assumed.

For a number `x`, `adapt_precision(T, x)` behaves as [`convert_real_type(T, x)`](@ref
convert_real_type) and `adapt_precision(T, typeof(x))` may be used to infer the
corresponding type with precision `T`.

Example:

```jldoctest; setup=:(using TypeUtils)
julia> adapt_precision(Float32, (1, 3//4, 3.0 - 2.0im, ("hello", false, 0x07, 1.0, π)))
(1, 0.75f0, 3.0f0 - 2.0f0im, ("hello", false, 0x07, 1.0f0, 3.1415927f0))
```

As can be seen, only non-integer numerical values are converted.

Basically, `adapt_precision` supports numbers, and arrays or tuples of numbers. It can be
specialized for other object types defined in foreign packages by specializing:

```julia
TypeUtils.adapt_precision(::Type{T}, x::SomeType) where {T<:TypeUtils.Precision} = ...
```

where `SomeType` is the object type and where the restriction `T<:TypeUtils.Precision` is to
make sure the above method is only called with a concrete floating-point type `T`.

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
# `convert_real_type` but leaving integers unchanged.
adapt_precision(::Type{T}, x::T) where {T<:Precision} = x
adapt_precision(::Type{T}, x::Integer) where {T<:Precision} = x
adapt_precision(::Type{T}, x::BareNumber) where {T<:Precision} = convert_real_type(T, x)
adapt_precision(::Type{T}, ::Type{T}) where {T<:Precision} = T
adapt_precision(::Type{T}, ::Type{S}) where {T<:Precision,S<:Integer} = S
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

Return a callable object `f` such that `f(x)` is equivalent to `adapt_precision(T, x)`. If
`T` is `AbstractFloat`, the default floating-point type
[`TypeUtils.default_precision`](@ref) is assumed.

"""
adapt_precision(::Type{AbstractFloat}) = adapt_precision(default_precision)
adapt_precision(::Type{T}) where {T<:Precision} = Converter(adapt_precision, T)
adapt_precision(::Type{T}) where {T} = throw_not_precision(T)

"""
    adapt_multiplier_precision(α, x)
    adapt_multiplier_precision(α, typeof(x))
    adapt_multiplier_precision(α, eltype(x))
    adapt_multiplier_precision(T::Type{<:Number}, α)

Adapt the precision of the multiplier `α` to compute `α*x` with `x` a numerical array. The
returned number has the same precision as type `T` of the elements of array `x`. If `α` is a
static number (see [`Neutrals.is_static_number`](@ref)) or if `T` has no concrete
floating-point precision, `α` is returned unchanged.

"""
adapt_multiplier_precision(α::Number, x::AbstractArray) =
    adapt_multiplier_precision(α, typeof(x))
adapt_multiplier_precision(α::Number, ::Type{T}) where {T<:AbstractArray} =
    adapt_multiplier_precision(α, eltype(T))
adapt_multiplier_precision(α::Number, ::Type{T}) where {T<:Number} =
    adapt_multiplier_precision(T, α)

adapt_multiplier_precision(::Type{T}, α::Number) where {T<:Number} =
    adapt_multiplier_precision(get_precision(T)::AbstractFloat, α)

adapt_multiplier_precision(::Type{T}, α::Number) where {T<:AbstractFloat} =
    (is_static_number(α) || !isconcretetype(T)) ? α : convert_floating_point_type(T, α)
