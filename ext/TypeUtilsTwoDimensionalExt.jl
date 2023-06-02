module TypeUtilsTwoDimensionalExt

if isdefined(Base, :get_extension)
    using TwoDimensional
    using TypeUtils
else
    using ..TwoDimensional
    using ..TypeUtils
end

for (X,N,T) in ((:(TwoDimensional.Point), 2, Real),
                (:(TwoDimensional.WeightedPoint), 3, AbstractFloat),
                (:(TwoDimensional.BoundingBox), 4, Real),)
    @eval begin
        TypeUtils.as(::Type{$X}, x::NTuple{$N,$T}) = $X(x...)
        TypeUtils.as(::Type{$X{T}}, x::NTuple{$N,$T}) where {T} = $X{T}(x...)
        TypeUtils.as(::Type{Tuple}, x::$X) = Tuple(x)
        TypeUtils.as(::Type{NTuple{$N,T}}, x::$X) where {T<:$T} = map(as(T), Tuple(x))
    end
end

end # module
