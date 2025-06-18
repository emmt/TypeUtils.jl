# Methods related to functions.

"""
    return_type(f, argtypes...) -> T

yields the type of the result returned by the callable object `f` when called
with arguments of types `argtypes...`.

See the warning in the documentation of `Base.promote_op` for the fragility of such
inference in some cases. There are no such issues if `f` is an instance of
[`TypeStableFunction`](@ref), e.g. built by [`as_return`][@ref), however `argtypes...` are
not checked for validity for such objects.

"""
return_type(::TypeStableFunction{T}, ::DataType...) where {T} = T
return_type(::Type{<:TypeStableFunction{T}}, ::DataType...) where {T} = T
return_type(f, argtypes::DataType...) = Base.promote_op(f, argtypes...)

"""
    as_return(T, f) -> g
    TypeStableFunction{T}(f) -> g
    TypeStableFunction(f, argtypes...) -> g

yield a callable object `g` that wraps callable `f` for guaranteed returned type `T`.
Alternatively, the type(s) `argtypes...` of the function argument(s) can be specified to
infer the returned type `T`. Then the following holds:

    g(args...; kwds...) === as(T, f(args...; kwds...))

for any arguments `args...` and keywords `kwds...`. Note that due to limitation of the
`Base.promote_op` method, it is currently not possible to infer `T` based on the types of
the keywords.

Methods [`return_type(g)`](@ref) and `parent(g)` may be used to retrieve `T` and `f`
respectively.

If `T` is specified, a similar object is given by:

    g = as(T)∘f

"""
as_return(::Type{T}, f) where {T} = TypeStableFunction{T}(f)

@doc @doc(as_return) TypeStableFunction

# Outer constructor.
function TypeStableFunction(f, argtypes::DataType...)
    T = Base.promote_op(f, argtypes...)
    R = T isa Union ? promote_type(ntuple(i -> getfield(T, i), Val(nfields(T)))...) : T
    isconcretetype(R) || throw(ArgumentError("cannot promote `$T` to a single concrete type"))
    return TypeStableFunction{R}(f)
end

# Conversion constructors.
TypeStableFunction{T}(f::TypeStableFunction{T}) where {T} = f
TypeStableFunction{T}(f::TypeStableFunction) where {T} = TypeStableFunction{T}(parent(f))

# Abstract constructor.
AbstractTypeStableFunction(f, argtypes::DataType...) = TypeStableFunction(f, argtypes...)
AbstractTypeStableFunction{T}(f) where {T} = TypeStableFunction{T}(f)

# Make instances of TypeStableFunction callable.
@inline (obj::TypeStableFunction{T})(args...; kwds...) where {T} =
    as(T, parent(obj)(args...; kwds...))

# Extend base methods.
Base.parent(obj::TypeStableFunction) = getfield(obj, :callable)

Base.return_types(::AbstractTypeStableFunction{T}; kwds...) where {T} = (T,)
Base.return_types(::AbstractTypeStableFunction{T}, ::DataType; kwds...) where {T} = (T,)

Base.promote_op(::AbstractTypeStableFunction{T}, ::DataType...) where {T} = T

for cls in (:AbstractTypeStableFunction, :TypeStableFunction,)
    @eval begin
        Base.convert(::Type{T}, f::T) where {T<:$(cls)} = f
        Base.convert(::Type{$(cls){T}}, f) where {T} = $(cls){T}(f)
    end
end
