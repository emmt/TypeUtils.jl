"""
    destructure(obj) -> vals::Tuple
    destructure(Tuple, obj) -> vals::Tuple
    destructure(Vector, obj) -> vals::Vector
    destructure(Vector{T}, obj) -> vals::Vector{T}

destructure object `obj` as a tuple or a vector of field values. Any structures in `obj`
are recursively destructured. For example, a complex is destructured as two reals.

See also [`destructure!`](@ref), [`restructure`](@ref), and [`struct_length`](@ref).

"""
destructure(::Type{Vector{T}}, obj) where {T} =
    destructure!(Vector{T}(undef, struct_length(obj)), obj)

destructure(::Type{Vector}, obj) = collect(destructure(obj))

destructure(::Type{Tuple}, obj) = destructure(obj)

@generated destructure(obj::T) where {T} = encode(destructure, :obj, T)

function encode(::typeof(destructure), base::Union{Symbol,Expr,QuoteNode}, ::Type{T}) where {T}
    expr = Expr(:tuple) # start with empty tuple
    encode!(destructure, expr.args, base, T)
    return quote
        $(Expr(:meta, :inline))
        return $expr
    end
end

function encode!(::typeof(destructure), code::AbstractVector,
                 base::Union{Symbol,Expr,QuoteNode}, ::Type{T}) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        for k in 1:fieldcount(T)
            encode!(destructure, code, :(getfield($base, $k)), fieldtype(T, k))
        end
    else
        push!(code, base)
    end
    nothing
end

"""
    destructure!(vals, obj; offset = firstindex(vals) - 1) -> vals

destructures object `obj` into `vals[offset+1:offset+n]` and returns `vals`. Here `n =
struct_length(obj)` is the total number of values stored by `obj`.

See also [`destructure`](@ref), [`restructure`](@ref), and
[`struct_length`](@ref).

"""
@generated function destructure!(vals::AbstractVector, obj::T;
                                 offset::Integer = firstindex(vals) - 1) where {T}
    return encode(destructure!, :vals, :offset, :obj, T)
end

function encode(::typeof(destructure!), arr::Symbol, off::Symbol,
                base::Union{Symbol,Expr,QuoteNode}, ::Type{T}) where {T}
    code = Expr(:block, Expr(:meta, :inline))
    encode!(destructure!, code.args,
            (i, x) -> :($arr[$off + $i] = $x), # function to store each field
            base, T, 0)
    push!(code.args, :(return $arr))
    return code
end

function encode!(::typeof(destructure!), code::AbstractVector, f,
                 base::Union{Symbol,Expr,QuoteNode}, ::Type{T}, i::Int) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        for k in 1:fieldcount(T)
            i = encode!(destructure!, code, f, :(getfield($base, $k)), fieldtype(T, k), i)
        end
    else
        i += 1
        push!(code, f(i, base))
    end
    return i
end

"""
    restructure(T, vals; offset = firstindex(vals) - 1) -> obj::T

restructures values `vals[offset+1:offset+n]` into an object `obj` of type `T`. Here `n =
struct_length(T)` is the total number of values stored by an object of type `T`.

The default constructor must exist for type `T` and, recursively, for any structured
fields of `T`.

See also [`destructure`](@ref), [`destructure!`](@ref), and [`struct_length`](@ref).

For an immutable concrete object `obj`, the following identity holds:

    restructure(typeof(obj), destructure(obj)) === obj

"""
@generated function restructure(::Type{T}, vals::Union{Tuple,AbstractVector};
                                offset::Integer = firstindex(vals) - 1) where {T}
    return encode(restructure, T, :vals, :offset)
end

function encode(::typeof(restructure), ::Type{T}, vals::Symbol, off::Symbol) where {T}
    code = Expr(:block, Expr(:meta, :inline))
    encode!(restructure, code, T, i -> :($vals[$off + $i]), 0)
    return code
end

function encode!(::typeof(restructure), code::Expr, ::Type{T}, f, i::Int) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        expr = if T <: Tuple
            Expr(:tuple)
        elseif isconcretetype(T)
            Expr(:call, T)
        else
            Expr(:call, parameterless(T))
        end
        for k in 1:fieldcount(T)
            i = encode!(restructure, expr, fieldtype(T, k), f, i)
        end
        if T <: Tuple
            push!(code.args, :(convert($T, $expr)))
        else
            push!(code.args, expr)
        end
    else
        i += 1
        push!(code.args, :($(f(i))))
    end
    return i
end

"""
    struct_length(x) -> n
    struct_length(typeof(x)) -> n

yield the total number of values stored by the fields of a structured object `x`. The
result only depends on the type of `x`.

"""
struct_length(x) = struct_length(typeof(x))

@generated function struct_length(::Type{T}) where {T}
    if isstructtype(T) && fieldcount(T) > 0
        n = 0
        for k in 1:fieldcount(T)
            n += struct_length(fieldtype(T, k))
        end
        return n
    else
        return 1
    end
end
