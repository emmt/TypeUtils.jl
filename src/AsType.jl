module AsType

export as

"""
    as(T, x)

yields `x` converted to type `T`.

"""
as(::Type{T}, x::T) where {T} = x
as(::Type{T}, x) where {T} = convert(T, x)::T

"""
    as(T)

yields a callable object which converts its argument to type `T`.

"""
as(::Type{T}) where {T} = As{T}()
struct As{T} <: Function; end
(::As{T})(x) where {T} = as(T, x)

end
