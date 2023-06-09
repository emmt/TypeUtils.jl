module TypeUtils

export
    as,
    as_eltype,
    as_return,
    convert_eltype,
    parameterless,
    promote_eltype,
    return_type

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
# https://stackoverflow.com/questions/42229901/getting-the-parameter-less-type
#
# NOTE: In old versions of Julia, the field name was `:primary`, but since
#       Julia 0.7, it should be `:wrapper`.
@inline parameterless(::Type{T}) where {T} = getfield(Base.typename(T), :wrapper)

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

yields an array `B` with the same entries as `A` except that their type is `T`.
If `T` is the element type of `A`, then `A` is returned.

!!! warning
    Calling this method for ranges yields a vector except if `T` is the element
    type of `A`. This is necessary to insure that `B` and `A` have the same
    size and that `B[i] == convert(T,A[i])` holds for all indices `i` in `A`.

"""
convert_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractArray) where {T} = convert(AbstractArray{T}, A)

"""
    as_eltype(T, A) -> B

yields an array which lazily converts its entries to type `T`. More
specifically, a call like `B[i]` yields `as(T,A[i])`.

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

function __init__()
    @static if !isdefined(Base, :get_extension)
        @require TwoDimensional="1907e7ba-7586-4310-a2ba-dd01462aeb50" include(
            "../ext/TypeUtilsTwoDimensionalExt.jl")
    end
end

end # module TypeUtils
