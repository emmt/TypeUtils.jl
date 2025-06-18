# Methods related to arrays.

"""
    as_array_axes(args...) -> rngs::ArrayAxes

converts array dimensions or ranges `args...` to a canonical form of array axes, that is a
tuple of `AbstractUnitRange{eltype(Dims)}`s. Any integer in `args...` is replaced by an
instance of `Base.OneTo{eltype(Dims)}`.

The array dimensions or ranges may also be provided as a tuple.
[`RelaxedArrayShape{N}`](@ref) is the union of types of `N`-tuples to which
`as_array_axes` is applicable.

Also see [`as_array_shape`](@ref), [`as_array_size`](@ref), [`as_array_axis`](@ref),
[`ArrayAxes`](@ref), `Dims`, and [`new_array`](@ref).

"""
as_array_axes(::Tuple{}) = ()
as_array_axes(args::eltype(RelaxedArrayShape)...) = as_array_axes(args)
as_array_axes(rngs::ArrayAxes) = rngs
as_array_axes(args::RelaxedArrayShape) = map(as_array_axis, args)

"""
    as_array_axis(arg) -> rng::ArrayAxis

converts array dimension or range `arg` to a canonical array axis, that is an instance of
`AbstractUnitRange{eltype(Dims)}`. If `arg` is an integer, `Base.OneTo{eltype(Dims)}(arg)`
is returned. [`eltype(RelaxedArrayShape)`](@ref RelaxedArrayShape) is the union of types
to which `as_array_axis` is applicable.

Also see [`as_array_axes`](@ref), [`as_array_dim`](@ref), and `Dims`.

"""
as_array_axis(dim::Integer) = Base.OneTo{Dim}(dim)
as_array_axis(rng::ArrayAxis) = rng
as_array_axis(rng::AbstractUnitRange{<:Integer}) = ArrayAxis(rng)
as_array_axis(rng::AbstractRange{<:Integer}) =
     isone(step(rng)) ? UnitRange{Dim}(first(rng), last(rng)) : throw_non_unit_step(rng)

"""
    as_array_dim(arg) -> dim::eltype(Dims)

converts array dimension or range `arg` to a canonical array dimension, that is an
`eltype(Dims)`. If `arg` is a unit-step range, its length is returned.
[`eltype(RelaxedArrayShape)`](@ref RelaxedArrayShape) is the union of types to which
`as_array_dim` is applicable.

Also see [`as_array_size`](@ref), [`as_array_axis`](@ref), and `Dims`.

"""
as_array_dim(dim::Dim) = dim
as_array_dim(dim::Integer) = as(Dim, dim)
as_array_dim(rng::AbstractUnitRange{<:Integer}) = as_array_dim(length(rng))
as_array_dim(rng::AbstractRange{<:Integer}) =
     isone(step(rng)) ? as_array_dim(length(rng)) : throw_non_unit_step(rng)

@noinline throw_non_unit_step(rng::AbstractRange) = throw(ArgumentError(
    "invalid non-unit step ($(step(rng))) for array axis"))

"""
    as_array_shape(args...) -> r::Union{Dims,ArraysAxes}

converts array dimensions or ranges `args...` to a canonical form of array shape, one of:

* array size, that is a tuple of `Int`s. This is the result if all of `args...` are
  integers or instances of `Base.OneTo`, the latter, if any, being replaced by their
  lengths.

* array axes, that is a tuple of `AbstractUnitRange{Int}`s. This is the result if any of
  `args...` are non-`Base.OneTo` ranges, the integers being converted to instances of
  `Base.OneTo{eltype(Dims)}`.

The array dimensions or ranges may also be provided as a tuple.
[`RelaxedArrayShape{N}`](@ref) is the union of types of `N`-tuples to which
`as_array_shape` is applicable.

Also see [`as_array_size`](@ref), [`as_array_axes`](@ref), [`ArrayAxes`](@ref), `Dims`,
and [`new_array`](@ref).

"""
as_array_shape(::Tuple{}) = ()
as_array_shape(args::eltype(RelaxedArrayShape)...) = as_array_shape(args)
as_array_shape(args::RegularArrayShape) = as_array_size(args)
as_array_shape(args::RelaxedArrayShape) = as_array_axes(args)

"""
    as_array_size(args...) -> dims::Dims

converts array dimensions or ranges `args...` to a canonical form of array size, that is a
tuple of `eltype(Dims)`s. Any range in `args...` is replaced by its length.

The array dimensions or ranges may also be provided as a tuple.
[`RelaxedArrayShape{N}`](@ref) is the union of types of `N`-tuples to which
`as_array_size` is applicable.

Also see [`as_array_shape`](@ref), [`as_array_axes`](@ref), [`as_array_dim`](@ref),
`Dims`, and [`new_array`](@ref).

"""
as_array_size(::Tuple{}) = ()
as_array_size(args::eltype(RelaxedArrayShape)...) = as_array_size(args)
as_array_size(dims::Dims) = dims
as_array_size(args::RelaxedArrayShape) = map(as_array_dim, args)

"""
    new_array(T, inds...) -> A
    new_array(T, (inds...,)) -> A

create a new array with element type `T` and shape defined by `inds...`, a list of array
dimension lengths and/or index ranges. The shape may also be specified as a tuple.

If `inds...` contains any index range other than `Base.OneTo`, an `OffsetArray{T}` is
returned; otherwise an `Array{T}` is returned. In the former case, an exception is thrown
if the package `OffsetArrays` has not been loaded.

Also see [`as_array_shape`](@ref), [`as_array_axes`](@ref), and [`as_array_size`](@ref).

"""
new_array(::Type{T}, shape::eltype(RelaxedArrayShape)...) where {T} = new_array(T, shape)
new_array(::Type{T}, shape::RelaxedArrayShape) where {T} = new_array(T, as_array_shape(shape))
new_array(::Type{T}, shape::RegularArrayShape) where {T} = new_array(T, as_array_size(shape))
new_array(::Type{T}, shape::Tuple{}) where {T} = Array{T}(undef)
new_array(::Type{T}, shape::Dims{N}) where {T,N} = Array{T,N}(undef, shape)
new_array(::Type,    shape::Union{Unsupported,ArrayAxes}) =
    error("package `OffsetArrays` must be loaded for such array index ranges")

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

converts the element type of object/type `A` to type `T`. The returned object/type is
similar to `A` except maybe for the element type. For example, if `A` is a range, then `B`
is also a range. If `T` is the element type of `A`, then `A` may be returned.

Consider using [`as_eltype(T, A)`](@ref as_eltype) to build an object that lazily performs
the conversion.

To simplify extending `convert_eltype` for objects `A` of given type, the default behavior
is:

    convert_eltype(T, A) = as(convert_eltype(T, typeof(A)), A)

so that it may be sufficient to extend `convert_eltype` for the type of the objects.

"""
convert_eltype(::Type{T}, x::X) where {T,X} = as(convert_eltype(T, X), x)
convert_eltype(::Type{T}, ::Type{X}) where {T,X} =
    error("don't know how to convert the element type of type `$X` to `$T`")

# As of Julia 1.11 `AbstractMatrix{T}(A)` or `AbstractArray{T}(A)` can be used to convert
# element type of `A` for Adjoint, Bidiagonal, Diagonal, Hermitian,
# LinearAlgebra.LQPackedQ, LowerTriangular, SymTridiagonal, Symmetric, Transpose,
# Tridiagonal, UnitLowerTriangular, UnitUpperTriangular, UpperHessenberg, UpperTriangular,
# etc.
convert_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
convert_eltype(::Type{T}, A::AbstractArray) where {T} = AbstractArray{T}(A)
convert_eltype(::Type{T}, ::Type{<:Array{<:Any,N}}) where {T,N} = Array{T,N}
convert_eltype(::Type{T}, ::Type{<:AbstractArray{<:Any,N}}) where {T,N} = AbstractArray{T,N}

# `LinearAlgebra.Factorization{T}(A)` can be used to convert element type of `A` for QR,
# LinearAlgebra.QRCompactWY, QRPivoted, LQ, Cholesky, CholeskyPivoted, LU, LDLt,
# BunchKaufman, SVD, etc.
convert_eltype(::Type{T}, A::Factorization{T}) where {T} = A
convert_eltype(::Type{T}, A::Factorization) where {T} = Factorization{T}(A)
if VERSION < v"1.7.0-beta1"
    # For old Julia versions, the above is not sufficient for SVD.
    convert_eltype(::Type{T}, A::SVD{T}) where {T} = A
    convert_eltype(::Type{T}, A::SVD) where {T} =
        SVD(convert_eltype(T, A.U), convert_eltype(real(T), A.S), convert_eltype(T, A.Vt))
end
convert_eltype(::Type{T}, A::Hessenberg{T}) where {T} = A
convert_eltype(::Type{T}, A::Hessenberg) where {T} = throw(ArgumentError(
    "changing element type of Hessenberg decomposition is not supported, consider re-computing the decomposition"))

# For `Adjoint` and `Transpose`, we want to preserve this structure.
for S in (:Adjoint, :Transpose)
    @eval begin
        convert_eltype(::Type{T}, A::$S{T}) where {T} = A
        convert_eltype(::Type{T}, A::$S) where {T} = $S(convert_eltype(T, parent(A)))
    end
end

# Convert element type for numbers.
convert_eltype(::Type{T}, ::Type{<:Number}) where {T} = T

# Convert element type for tuples. See `_countuple` in `base/tuple.jl` for the best
# way to extract the number of elements in a tuple given its type.
convert_eltype(::Type{T}, ::Type{<:NTuple{N,Any}}) where {N,T} = NTuple{N,T}
convert_eltype(::Type{T}, A::NTuple{N,T}) where {N,T} = A
convert_eltype(::Type{T}, A::Tuple) where {T} = map(as(T), A)

# Convert element type for `Base.OneTo` and `UnitRange`. For `T` non-integer,
# a `Base.OneTo` instance becomes a `UnitRange` one (in a predictible way).
convert_eltype(::Type{T}, ::Type{<:OneTo}) where {T<:Integer} = OneTo{T}
convert_eltype(::Type{T}, ::Type{<:OneTo}) where {T} = UnitRange{T}
convert_eltype(::Type{T}, ::Type{<:UnitRange}) where {T} = UnitRange{T}
convert_eltype(::Type{T}, r::Union{OneTo,UnitRange}) where {T} =
    as(convert_eltype(T, typeof(r)), r)

# Convert element type for AbstractUnitRange{T} <: OrdinalRange{T,T}.
convert_eltype(::Type{T}, ::Type{<:AbstractUnitRange}) where {T} = typeof(zero(T):oneunit(T))
convert_eltype(::Type{T}, r::AbstractUnitRange{T}) where {T} = r
convert_eltype(::Type{T}, r::AbstractUnitRange) where {T} = as(T, first(r)):as(T, last(r))

# Convert element type for other range types.
convert_eltype(::Type{T}, ::Type{<:AbstractRange}) where {T} = AbstractRange{T}
convert_eltype(::Type{T}, r::AbstractRange{T}) where {T} = r
convert_eltype(::Type{T}, r::AbstractRange) where {T} =
    as(T, first(r)):as(T, step(r)):as(T, last(r))

# Convert element type for LinRange{T,L<:Integer} <: AbstractRange{T}.
convert_eltype(::Type{T}, r::LinRange{T}) where {T} = r
convert_eltype(::Type{T}, r::LinRange) where {T} =
    LinRange(as(T, first(r)), as(T, last(r)), length(r))

"""
    convert_eltype(T) -> f

yields a callable object `f` such that `f(x)` yields `convert_eltype(T, x)` for any `x`.

"""
convert_eltype(::Type{T}) where {T} = Converter(convert_eltype, T)

"""
    as_eltype(T, A) -> B

yields an array which lazily converts its entries to type `T`. More specifically, a call
like `B[i]` yields `as(T,A[i])`.

Consider using [`convert_eltype(T, A)`](@ref convert_eltype) to perform the conversion
once and immediately.

"""
as_eltype(::Type{T}, A::AbstractArray{T}) where {T} = A
as_eltype(::Type{T}, A::AbstractArray) where {T} = AsEltype{T}(A)

Base.parent(A::AsEltype) = getfield(A, :parent)

# Implement abstract array API for `AsEltype` objects.
for func in (:axes, :length, :size)
    @eval Base.$func(A::AsEltype) = $func(parent(A))
end
for (L, S, Idecl, Icall) in ((false, :IndexCartesian, :(I::Vararg{Int,N}), :(I...)),
                             (true,  :IndexLinear,    :(i::Int),           :(i)))
    @eval begin
        Base.IndexStyle(::Type{<:AsEltype{T,N,$L}}) where {T,N} = $S()
        @inline function Base.getindex(A::AsEltype{T,N,$L}, $Idecl) where {T,N}
            @boundscheck checkbounds(A, $Icall)
            r = @inbounds getindex(parent(A), $Icall)
            return as(T, r)
        end
        @inline function Base.setindex!(A::AsEltype{T,N,$L}, x, $Idecl) where {T,N}
            @boundscheck checkbounds(A, $Icall)
            @inbounds setindex!(parent(A), x, $Icall)
            return A
        end
    end
end

Base.similar(A::AsEltype, ::Type{T}) where {T} = similar(parent(A), T)
for shape in (:Dims,
              :(Tuple{Integer,Vararg{Integer}}),
              :(Tuple{Union{Integer,UnitRange{<:Integer}},
                      Vararg{Union{Integer,UnitRange{<:Integer}}}}))
    @eval Base.similar(A::AsEltype, ::Type{T}, shape::$shape) where {T} =
        similar(parent(A), T, shape)
end
