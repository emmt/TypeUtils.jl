module TypeUtilsUnitfulExt
if isdefined(Base, :get_extension)
    using TypeUtils, Unitful
    using Unitful: AbstractQuantity, Quantity, ustrip, unit
else
    using ..TypeUtils, ..Unitful
    using ..Unitful: AbstractQuantity, Quantity, ustrip, unit
end


# Extend bare_type, real_type, convert_bare_type, and convert_real_type.
for (f, g, S) in ((:(TypeUtils.bare_type), :(TypeUtils.convert_bare_type), :(TypeUtils.BareNumber)),
                  (:(TypeUtils.real_type), :(TypeUtils.convert_real_type), :Real))
    @eval begin
        $f(::Type{<:AbstractQuantity{T}}) where {T} = $f(T)
        $g(::Type{T}, x::AbstractQuantity{T}) where {T<:$S} = x
        @inline $g(::Type{T}, x::AbstractQuantity) where {T<:$S} =
            $g(T, ustrip(x))*unit(x)
        $g(::Type{T}, ::Type{Quantity{T,D,U}}) where {T<:$S,D,U} =
            Quantity{T,D,U}
        @inline $g(::Type{T}, ::Type{Quantity{R,D,U}}) where {T<:$S,R,D,U} =
            Quantity{$g(T,R),D,U}
    end
end

# Extend unitless (only needed for values).
TypeUtils.unitless(x::AbstractQuantity) = ustrip(x)

end # module
