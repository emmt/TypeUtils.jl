module AsType

export as

using Requires

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

function __init__()
    @require TwoDimensional="1907e7ba-7586-4310-a2ba-dd01462aeb50" begin
        for (X,N,T) in ((:(TwoDimensional.Point), 2, Real),
                        (:(TwoDimensional.WeightedPoint), 3, AbstractFloat),
                        (:(TwoDimensional.BoundingBox), 4, Real),)
            @eval begin
                as(::Type{$X}, x::NTuple{$N,$T}) = $X(x...)
                as(::Type{$X{T}}, x::NTuple{$N,$T}) where {T} = $X{T}(x...)
                as(::Type{Tuple}, x::$X) = Tuple(x)
                as(::Type{NTuple{$N,T}}, x::$X) where {T<:$T} = map(as(T), Tuple(x))
            end
        end
    end
end

end
