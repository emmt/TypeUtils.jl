"""

`TypeUtils.LazyMaps` provides a single method, `lazymap`, to build lazily mapped arrays or
collections for Julia. In a nutshell:

```julia
B = lazymap([T::Type,] f, A)
B = lazymap([T::Type,] f, A::AbstractArray, f_inv = inverse(f))
```

yields a view of `A` such that the `i`-th element of `B` is given by `Bᵢ = as(T, f(Aᵢ))`
with `Aᵢ` the `i`-th element of `A`. If `T` is unspecified, it is inferred from `f` and from
the element type of `A`. If `A` is an array, then `f_inv` is assumed to be the inverse of
`f` and `B[i] = x` amounts to do `A[i] = f_inv(x)`.

"""
module LazyMaps

export
    LazyMap,
    lazymap

using ..TypeUtils
using ..TypeUtils: @public
@public result

using InverseFunctions

struct LazyMapArray{T,N,F,A<:AbstractArray,L,I} <: AbstractArray{T,N}
    f::F     # callable
    arg::A   # input array
    f_inv::I # inverse function (callable or unknown)
    LazyMapArray{T}(f::F, arg::A, f_inv::I) where {T,N,F,A<:AbstractArray{<:Any,N},I} =
        new{T,N,F,A,IndexStyle(A) isa IndexLinear,I}(f, arg, f_inv)
end

struct LazyMapAny{T,N,F,A}
    f::F   # callable
    arg::A # input collection
    LazyMapAny{T}(f::F, arg::A) where {T,F,A} =
        new{T,infer_ndims(Base.IteratorSize(A)),F,A}(f, arg)
end

"""
    LazyMap{T,N,F,A}

Union of lazy mapped arrays or collections. `T` is the element type, `N` the number of
dimensions, `F` the function type, and `A` the type of the argument of the lazy map.

"""
const LazyMap{T,N,F,A} = Union{LazyMapArray{T,N,F,A},
                               LazyMapAny{T,N,F,A}}

# Singleton type to indicate unknown parameter or type.
struct Unknown end

"""
    B = lazymap([T::Type,] f, A)
    B = lazymap([T::Type,] f, A::AbstractArray, f_inv = inverse(f))

Build a view of the array or iterator `A` such that the `i`-th element of `B` is given by
`Bᵢ = as(T, f(Aᵢ))` with `Aᵢ` the `i`-th element of `A`.

Optional argument `T` is to explicitly specify the element type of `B`; otherwise, it is
inferred from `f` and from the element type of `A`. The lazy map `B` has type-stable element
type in the sense that its element have guaranteed type `T`, even though `T` may be
abstract.

If `A` is an array, `f_inv` is the assumed inverse of `f` such that `B[i] = x` has the side
effect of modifying `A` by `A[i] = f_inv(x)`. If unspecified, `f_inv` is inferred by the
`inverse` method of the
[`InverseFunctions`](https://github.com/JuliaMath/InverseFunctions.jl) package. If `f_inv =
throw`, a read-only lazy map array is returned even though `inverse(f)` is known. Similarly,
if `f = throw`, a write-only lazy map object is returned (you probably want to specify
`f_inv` in this case).

As a special case:

    C = lazymap(T::Type, A)

builds an object `C` that lazily maps the **constructor** `T` to the elements of `A`. This
is not exactly the same as:

    B = lazymap(T::Type, identity, A)

which builds an object `B` that lazily **converts** the elements of `A` to type `T`. In
other words, the `i`-th element of `C` is given by `Cᵢ = T(Aᵢ)::T`, while the `i`-th element
of `B` is given by `Bᵢ = as(T, Aᵢ)`. In both cases, it is asserted that `Cᵢ` and `Bᵢ` are of
type `T`. The two are equivalent if `T` is a numeric type (a sub-type of `Number`).

## See also

The [`as`](@ref TypeUtils.as) and `inverse` methods respectively provided by the
 [`TypeUtils`](https://github.com/emmt/TypeUtils.jl) and
 [`InverseFunctions`](https://github.com/JuliaMath/InverseFunctions.jl) packages.

"""
lazymap(f, arg::Any) = lazymap(infer_eltype(f, arg), f, arg)
lazymap(f, arg::Any, f_inv) = lazymap(infer_eltype(f, arg), f, arg, f_inv)

lazymap(::Type{T}, arg::AbstractArray) where {T} = lazymap(T, pass, arg, pass)
lazymap(::Type{T}, arg::Any) where {T} = lazymap(T, pass, arg)
lazymap(::Type{T}, f, arg::AbstractArray) where {T} = lazymap(T, f, arg, inverse(f))
lazymap(::Type{T}, f, arg::AbstractArray, f_inv) where {T} = LazyMapArray{T}(f, arg, f_inv)
lazymap(::Type{T}, f, arg::Any) where {T} = LazyMapAny{T}(f, arg)
@noinline lazymap(::Type{T}, f, arg::A, f_inv) where {T,A} = throw(ArgumentError(
    "in `lazymap([T::Type,] f, arg, f_inv)`, `arg` must be an array, got an instance of `$A`"))

infer_eltype(f, arg::Any) = _infer_eltype(Base.IteratorEltype(arg), f, arg)
_infer_eltype(trait::Base.HasEltype, f::typeof(throw), arg) = Unknown
_infer_eltype(trait::Base.HasEltype, f, arg) = Base.promote_op(f, eltype(arg))
_infer_eltype(trait::Base.IteratorEltype, f, arg) = Unknown

# For collections, the shape traits are inferred according to the rules for tuples: if
# `IteratorSize(A)` yields `HasShape{N}()`, then `A` has a length, a number of dimensions,
# a size, and axes. Otherwise, if `IteratorSize(A)` yields `HasLength()`, then `A` has a
# length but no number of dimensions, size, nor axes.
infer_ndims(trait::Base.HasShape{N}) where {N} = N
infer_ndims(trait::Base.IteratorSize) = Unknown

# Dummy function for lazy maps `B = lazymap(T, A)` computing their output as `T(x)::T`, not
# as `as(T, B.f(x))` which, in most cases, is similar to `convert(T, B.f(x))`. This function
# behaves like `identity` but has its own type. Using it to implement the `T(x)::T` behavior
# results in a smaller size for `B` which only stores one reference (to the collection
# argument) instead of 2 (to the collection and to `T`).
pass(x) = x
InverseFunctions.inverse(::typeof(pass)) = pass

# Accessors.
Base.parent(A::LazyMap) = getfield(A, :arg)
func(A::LazyMap) = getfield(A, :f)
inv_func(A::LazyMapArray) = getfield(A, :f_inv)

# Abstract array API for instances of LazyMapArray.

Base.length(m::LazyMapArray) = length(parent(m))
Base.size(m::LazyMapArray) = size(parent(m))
Base.axes(m::LazyMapArray) = axes(parent(m))
Base.strides(m::LazyMapArray) = strides(parent(m))
Base.stride(m::LazyMapArray, k::Integer) = stride(parent(m), k)

for shape in (:Dims,
              :(Tuple{Integer,Vararg{Integer}}),
              :(Tuple{Union{Integer,UnitRange{<:Integer}},
                      Vararg{Union{Integer,UnitRange{<:Integer}}}}))
    @eval Base.similar(m::LazyMapArray, ::Type{T}, shape::$shape) where {T} =
        similar(parent(m), T, shape)
end

for (style, (Idecl, Icall)) in (:IndexLinear    => (:(i::Int),           :(i)),
                                :IndexCartesian => (:(I::Vararg{Int,N}), :(I...)))
    linear = (style === :IndexLinear)
    @eval begin
        Base.IndexStyle(::Type{<:LazyMapArray{T,N,F,A,$linear}}) where {T,N,F,A} = $style()
        @inline function Base.getindex(m::LazyMapArray{T,N,F,A,$linear},
                                       $Idecl) where {T,N,F,A}
            @boundscheck checkbounds(m, $Icall)
            F === typeof(throw) && throw_write_only()
            x = @inbounds getindex(parent(m), $Icall)
            return result(m, x)
        end
        @inline function Base.setindex!(m::LazyMapArray{T,N,F,A,$linear,Finv}, x,
                                        $Idecl) where {T,N,F,A,Finv}
            @boundscheck checkbounds(parent(m), $Icall)
            Finv === typeof(throw) && throw_read_only()
            @inbounds setindex!(parent(m), inv_func(m)(x), $Icall)
            return m
        end
    end
end

"""
    TypeUtils.LazyMaps.result(B, Aᵢ) -> Bᵢ

Compute the value returned by the lazy map `B = lazymap([T,] f, A)`. Here, `Aᵢ` and `Bᵢ` are
the respective `i`-th value of `A` and `B`.

This method may be specialized based on the type of the callable `f` with the following
signature:

    TypeUtils.LazyMaps.result(B::LazyMap{T,N,typeof(f)}, Bᵢ) where {T,N} = ...

"""
result(m::LazyMap{T,N,F,A}, x) where {T,N,F,A} = as(T, func(m)(x))
result(m::LazyMap{T,N,typeof(pass),A}, x) where {T,N,A} = T(x)::T

# Iterator and (partial) abstract array API for instances of LazyMapAny.
Base.IteratorEltype(::Type{<:LazyMapAny{T,N,F,A}}) where {T,N,F,A} =
    Base.HasEltype()
Base.IteratorEltype(::Type{<:LazyMapAny{T,N,F,A}}) where {T<:Unknown,N,F,A} =
    Base.EltypeUnknown()

Base.IteratorSize(::Type{<:LazyMapAny{T,N,F,A}}) where {T,N,F,A} = Base.IteratorSize(A)

Base.eltype(m::LazyMapAny) = eltype(typeof(m))
Base.eltype(::Type{<:LazyMapAny{T}}) where {T} = T
Base.eltype(::Type{<:LazyMapAny{Unknown}}) = throw_unknown_eltype()

Base.ndims(m::LazyMapAny) = ndims(typeof(m))
Base.ndims(::Type{<:LazyMapAny{T,N}}) where {T,N} = N
Base.ndims(::Type{<:LazyMapAny{T,Unknown}}) where {T} = throw_unknown_ndims()

Base.length(m::LazyMapAny) = _length(Base.IteratorSize(m), parent(m))
_length(trait::Base.HasLength, arg) = length(arg)
_length(trait::Base.HasShape, arg) = prod(_size(trait, arg))
_length(trait::Base.IteratorSize, arg) = throw_unknown_length()

Base.size(m::LazyMapAny) = _size(Base.IteratorSize(m), parent(m))
_size(trait::Base.HasShape, arg) = map(length, axes(arg))
_size(trait::Base.IteratorSize, arg) = throw_unknown_shape()

Base.axes(m::LazyMapAny) = _axes(Base.IteratorSize(m), parent(m))
_axes(trait::Base.HasShape, arg) = axes(arg)
_axes(trait::Base.IteratorSize, arg) = throw_unknown_shape()

# Make an instance of LazyMapAny an iterable.
Base.iterate(m::LazyMapAny) = _iterate_result(m, iterate(parent(m)))
Base.iterate(m::LazyMapAny, s) = _iterate_result(m, iterate(parent(m), s))
_iterate_result(m::LazyMapAny{T}, r::Nothing) where {T} = nothing
_iterate_result(m::LazyMapAny{T}, r::Tuple{Any,Any}) where {T} =
    (result(m, r[1]), r[2])

@noinline throw_read_only() =
    throw(ArgumentError("attempt to write read-only lazily mapped object"))
@noinline throw_write_only() =
    throw(ArgumentError("attempt to read write-only lazily mapped object"))
@noinline throw_unknown_ndims() =
    throw(ArgumentError("collection in lazy map has no defined number of dimensions"))
@noinline throw_unknown_length() =
    throw(ArgumentError("collection in lazy map has no defined length"))
@noinline throw_unknown_shape() =
    throw(ArgumentError("collection in lazy map has no defined shape"))
@noinline throw_unknown_eltype() =
    throw(ArgumentError("collection in lazy map has no defined element type"))

end
