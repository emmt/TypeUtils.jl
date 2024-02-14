module TypeUtils

export
    as,
    as_eltype,
    as_return,
    convert_eltype,
    destructure,
    destructure!,
    parameterless,
    promote_eltype,
    restructure,
    return_type,
    struct_length

using Base: OneTo

if !isdefined(Base, :get_extension)
    using Requires
end

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
        # convert to the requested tuple type.
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
    f = as(T)

yields a callable object which converts its argument to type `T`. More
specifically, a call like `f(x)` yields `as(T, x)`.

"""
as(::Type{T}) where {T} = As{T}()

struct As{T} <: Function; end

(::As{T})(x) where {T} = as(T, x)

"""
    g = as_return(T, f)

yields a callable object such that `g(args...; kwds...)` returns `f(args...;
kwds...)` converted to type `T`. Methods [`return_type(g)`](@ref) and
`parent(g)` may be used to retrieve `T` and `f` respectively.

A similar object is given by:

    g = as(T)∘f

""" as_return

struct AsReturn{T,F}
    func::F

    # Inner contructor.
    AsReturn{T}(func::F) where {T,F} = new{T,F}(func)

    # Avoid multiple wrapping.
    AsReturn{T}(func::AsReturn{T}) where {T} = func
    AsReturn{T}(func::AsReturn) where {T} = AsReturn{T}(parent(func))
end

(obj::AsReturn{T})(args...; kwds...) where {T} = as(T, parent(obj)(args...; kwds...))

as_return(::Type{T}, func) where {T} = AsReturn{T}(func)

Base.parent(obj::AsReturn) = getfield(obj, :func)
Base.return_types(obj::AsReturn{T}) where {T} = (T,)
Base.promote_op(obj::AsReturn{T}, argtypes::Type...) where {T} = T

"""
    return_type(f, argtypes...) -> T

yields the type of the result returned by the callable object `f` when called
with arguments of types `argtypes...`.

See the warning in the documentation of `Base.promote_op` for the fragility of
such inference in some cases. There are no such issues if `f` is an object
built by [`as_return`][@ref), however `argtypes...` are not checked for
validity for such objects.

"""
return_type(obj::AsReturn{T}, argtypes::Type...) where {T} = T
return_type(::Type{<:AsReturn{T}}, argtypes::Type...) where {T} = T
return_type(f, argtypes::Type...) = Base.promote_op(f, argtypes...)

"""
    parameterless(T)

yields the type `T` without parameter specifications. For example:

```julia
julia> parameterless(Vector{Float32})
Array
```

""" parameterless
#
# This subject is discussed in different places:
#
# - https://discourse.julialang.org/t/typeutils-dealing-with-types-in-julia/101584/13
#
# - https://stackoverflow.com/questions/42229901/getting-the-parameter-less-type
#
# The latter leads to define the method as (in old versions of Julia, the field
# name was `:primary`, but since Julia 0.7, it should be `:wrapper`):
#
#     @inline parameterless(::Type{T}) where {T} = getfield(Base.typename(T), :wrapper)
#
# The actual implementation is borrowed from `ConstructionBase` and should be
# less likely to be broken by internal changes in Julia:
#
parameterless(::Type{T}) where {T} = getfield(parentmodule(T), nameof(T))

"""
    promote_eltype(args...)

yields the promoted element type of its arguments. Arguments `args...` may be
anything implementing the `eltype` method.

"""
promote_eltype() = promote_type()
promote_eltype(arg) = eltype(arg)
@inline promote_eltype(args...) = promote_type(map(eltype, args)...)

"""
    convert_eltype(T, A) -> B

converts the elements of `A` to type `T`. The returned object is similar to `A`
except maybe for the element type. For example, if `A` is a range, then `B` is
also a range. If `T` is the element type of `A`, then `A` is returned.

Consider using [`as_eltype(T, A)`](@ref) to build an object that lazily
performs the conversion.

"""
convert_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractArray) where {T} = convert(AbstractArray{T}, A)

# Convert element type for tuples.
convert_eltype(::Type{T}, A::NTuple{N,T}) where {N,T} = A
convert_eltype(::Type{T}, A::Tuple) where {T} = map(as(T), A)

# Convert element type for Base.OneTo{T<:Integer} <: AbstractUnitRange{T}.
# Conversion to non-integer element types is handled by the more general method
# for AbstractUnitRange.
convert_eltype(::Type{T}, A::OneTo{T}) where {T<:Integer} = A
convert_eltype(::Type{T}, A::OneTo) where {T<:Integer} = OneTo{T}(last(A))

# Convert element type for AbstractUnitRange{T} <: OrdinalRange{T,T}.
convert_eltype(::Type{T}, A::AbstractUnitRange{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractUnitRange) where {T} =
    as(T, first(A)):as(T, last(A))

# Convert element type for other range types.
convert_eltype(::Type{T}, A::AbstractRange{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractRange) where {T} =
    as(T, first(A)):as(T, step(A)):as(T, last(A))

# Convert element type for LinRange{T,L<:Integer} <: AbstractRange{T}.
convert_eltype(::Type{T}, A::LinRange{T}) where {T} = A
convert_eltype(::Type{T}, A::LinRange) where {T} =
    LinRange(as(T, first(A)), as(T, last(A)), length(A))

"""
    as_eltype(T, A) -> B

yields an array which lazily converts its entries to type `T`. More
specifically, a call like `B[i]` yields `as(T,A[i])`.

Consider using [`convert_eltype(T, A)`](@ref) to perform the conversion once
and immediately.

"""
as_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
as_eltype(::Type{T}, A::AbstractArray) where {T} = AsEltype{T}(A)

struct AsEltype{T,N,L,A<:AbstractArray} <: AbstractArray{T,N}
    parent::A
end

AsEltype{T}(arr::A) where {T,N,A<:AbstractArray{<:Any,N}} =
    AsEltype{T,N,IndexStyle(A)===IndexLinear(),A}(arr)

Base.parent(A::AsEltype) = A.parent

# Implement abstract array API for `AsEltype` objects.
for func in (:axes, :length, :size)
    @eval Base.$func(A::AsEltype) = $func(parent(A))
end
Base.IndexStyle(::Type{<:AsEltype{<:Any,<:Any,true}}) = IndexLinear()
Base.IndexStyle(::Type{<:AsEltype{<:Any,<:Any,false}}) = IndexCartesian()

@inline function Base.getindex(A::AsEltype{T,N,true}, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds r = getindex(parent(A), i)
    return as(T, r)
end

@inline function Base.getindex(A::AsEltype{T,N,false}, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds r = getindex(parent(A), I...)
    return as(T, r)
end

@inline function Base.setindex!(A::AsEltype{T,N,true}, x, i::Int) where {T,N}
    @boundscheck checkbounds(A, i)
    @inbounds setindex!(parent(A), x, i)
    return A
end

@inline function Base.setindex!(A::AsEltype{T,N,false}, x, I::Vararg{Int,N}) where {T,N}
    @boundscheck checkbounds(A, I...)
    @inbounds setindex!(parent(A), x, I...)
    return A
end

"""
    destructure(obj) -> vals::Tuple

destructures object `obj` as a tuple of field values. Any structures in `obj`
are recursively destructured.

See also [`destructure!`](@ref), [`restructure`](@ref), and
[`struct_length`](@ref).

"""
function destructure end

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

destructures object `obj` into `vals[offset+1:offset+n]` and returns `vals`. Here
`n = struct_length(obj)` is the total number of values stored by `obj`.

See also [`destructure`](@ref), [`restructure`](@ref), and
[`struct_length`](@ref).

"""
function destructure! end

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

restructures values `vals[offset+1:offset+n]` into an object `obj` of type `T`.
Here `n = struct_length(T)` is the total number of values stored by an object
of type `T`.

The default constructors must exist for `T` and, recursively, for any
structured fields of `T`.

See also [`destructure`](@ref), [`destructure!`](@ref), and
[`struct_length`](@ref).

For an immutable concrete object `obj`, the following identity holds:

    restructure(typeof(obj), destructure(obj)) === obj

"""
function restructure end

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
    struct_length(T) -> n

yields the total number of values stored by the fields of a structured object
of type `T`. Argument may also be an object of type `T`.

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

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require TwoDimensional="1907e7ba-7586-4310-a2ba-dd01462aeb50" include(
            "../ext/TypeUtilsTwoDimensionalExt.jl")
    end
end

end # module TypeUtils
