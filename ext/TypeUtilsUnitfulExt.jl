module TypeUtilsUnitfulExt
if isdefined(Base, :get_extension)
    using TypeUtils, Unitful
else
    using ..TypeUtils, ..Unitful
end

# Extend bare_type, real_type, convert_bare_type, and convert_real_type.
for (f, g, S) in ((:(TypeUtils.bare_type), :(TypeUtils.convert_bare_type), :(TypeUtils.BareNumber)),
                  (:(TypeUtils.real_type), :(TypeUtils.convert_real_type), :Real))
    @eval begin
        $f(::Type{<:Unitful.AbstractQuantity{T}}) where {T} = $f(T)
        $g(::Type{T}, x::Unitful.AbstractQuantity{T}) where {T<:$S} = x
        @inline $g(::Type{T}, x::Unitful.AbstractQuantity) where {T<:$S} =
            $g(T, Unitful.ustrip(x))*Unitful.unit(x)
        $g(::Type{T}, ::Type{Unitful.Quantity{T,D,U}}) where {T<:$S,D,U} =
            Unitful.Quantity{T,D,U}
        @inline $g(::Type{T}, ::Type{Unitful.Quantity{R,D,U}}) where {T<:$S,R,D,U} =
            Unitful.Quantity{$g(T,R),D,U}
    end
end

# Extend unitless (only needed for values).
TypeUtils.unitless(x::Unitful.AbstractQuantity) = Unitful.ustrip(x)

end # module
