module TypeUtils

export as

using Requires

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
        # convert to the specific tuple type.
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

end # module TypeUtils
