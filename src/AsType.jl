module AsType

export as

"""
    as(T, x)

yields `x` converted to type `T`.

"""
as(::Type{T}, x::T) where {T} = x
as(::Type{T}, x) where {T} = convert(T, x)::T

# Use the constructor to convert to Cartesian index.
as(::Type{T}, x::T) where {T<:CartesianIndex} = x
as(::Type{CartesianIndex}, x::Tuple{Vararg{Integer}}) = CartesianIndex(x)
as(::Type{CartesianIndex{N}}, x::NTuple{N,Integer}) where {N} = CartesianIndex{N}(x)

# Convert Cartesian indices to tuples.
as(::Type{Tuple}, x::CartesianIndex) = Tuple(x)
as(::Type{Tuple}, x::CartesianIndices) = x.indices

# Conversion between symbols and strings is not supported by `convert`.
as(::Type{String}, x::Symbol) = String(x)
as(::Type{Symbol}, x::String) = Symbol(x)

"""
    as(T)

yields a callable object which converts its argument to type `T`.

"""
as(::Type{T}) where {T} = As{T}()
struct As{T} <: Function; end
(::As{T})(x) where {T} = as(T, x)

end
