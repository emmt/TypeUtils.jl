module TypeUtilsUnitfulExt
if isdefined(Base, :get_extension)
    using TypeUtils, Unitful
    using TypeUtils: BareNumber
    using Unitful: AbstractQuantity, Quantity, ustrip, unit
else
    using ..TypeUtils, ..Unitful
    using ..TypeUtils: BareNumber
    using ..Unitful: AbstractQuantity, Quantity, ustrip, unit
end


# Extend bare_type, real_type, convert_bare_type, and convert_real_type.
for (f, g, S) in ((:(TypeUtils.bare_type), :(TypeUtils.convert_bare_type), :BareNumber),
                  (:(TypeUtils.real_type), :(TypeUtils.convert_real_type), :Real))
    @eval begin
        $f(::Type{<:AbstractQuantity{T}}) where {T} = $f(T)
        $g(::Type{T}, x::AbstractQuantity{T}) where {T<:$S} = x
        $g(::Type{T}, x::AbstractQuantity   ) where {T<:$S} = $g(T, ustrip(x))*unit(x)
        $g(::Type{T}, ::Type{Quantity{T,D,U}}) where {T<:$S,D,U} = Quantity{T,D,U}
        $g(::Type{T}, ::Type{Quantity{R,D,U}}) where {T<:$S,D,U,R} = Quantity{$g(T,R),D,U}
    end
end

# Extend unitless (only needed for values).
TypeUtils.unitless(x::AbstractQuantity) = ustrip(x)

TypeUtils.get_precision(::Type{<:AbstractQuantity{T}}) where {T} = get_precision(T)

TypeUtils.adapt_precision(::Type{T}, x::Quantity{T,D,U}) where {T<:Precision,D,U} = x
TypeUtils.adapt_precision(::Type{T}, x::Quantity{S,D,U}) where {T<:Precision,D,U,S} =
    Quantity{adapt_precision(T, S), D, U}(x)
TypeUtils.adapt_precision(::Type{T}, ::Type{Quantity{S,D,U}}) where {T<:Precision,S,D,U} =
    Quantity{adapt_precision(T, S), D, U}

end # module
