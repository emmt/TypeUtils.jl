module TypeUtilsOffsetArraysExt

if isdefined(Base, :get_extension)
    using TypeUtils, OffsetArrays
else
    using ..TypeUtils, ..OffsetArrays
end

# Extend `new_array` for offset axes.
TypeUtils.new_array(::Type{T}, rngs::ArrayAxes{N}) where {T,N} =
    OffsetArray(Array{T}(undef, as_array_size(rngs)), rngs)

end # module
