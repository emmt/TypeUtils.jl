# Dealing with types in Julia

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](./LICENSE.md)
[![Build Status](https://github.com/emmt/TypeUtils.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/TypeUtils.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/TypeUtils.jl?svg=true)](https://ci.appveyor.com/project/emmt/TypeUtils-jl)
[![version](https://juliahub.com/docs/General/TypeUtils/stable/version.svg)](https://juliahub.com/ui/Packages/General/TypeUtils)
[![Coverage](https://codecov.io/gh/emmt/TypeUtils.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/TypeUtils.jl)
[![deps](https://juliahub.com/docs/General/TypeUtils/stable/deps.svg)](https://juliahub.com/ui/Packages/General/TypeUtils?t=2)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Package `TypeUtils` provides useful methods to deal with types in
[Julia](https://www.julialang.org) and facilitate coding with numbers whether they have
units or not. The package provides methods to strip units from numbers or numeric types,
convert the numeric type of quantities (preserving their units if any), determine
appropriate numeric type to carry computations mixing numbers with different types and/or
units. These methods make it easy to write code that works consistently for numbers with any
units (including none). The intention is that the `TypeUtils` package automatically extends
its exported methods when packages such as
[`Unitful`](https://github.com/PainterQubits/Unitful.jl) are loaded.


## Cast value to type

The method, `as` is designed to *cast* a value to a given type. The name was inspired by
the built-in [Zig](https://ziglang.org/) function
[`@as`](https://ziglang.org/documentation/master/#as).

A first usage is:

``` julia
as(T, x)
```

which yields `x` converted to type `T`. This behaves like `convert(T,x)::T` doing nothing
if `x` is already of type `T` and performing the conversion and the type assertion
otherwise.

By default, the `as` method calls `convert` only if needed but also implements a number of
conversions not supported by `convert`. The `as` method is therefore a bit more versatile
than `convert` while relaxing the bother to remember which function or constructor to call
to efficiently perform the intended conversion. For example:

``` julia
julia> as(Tuple, CartesianIndex(1,2,3)) # yields tuple of indices
(1, 2, 3)

julia> as(CartesianIndex, (1,2,3)) # calls constructor
CartesianIndex(1, 2, 3)

julia> as(Tuple, CartesianIndices(((-2:5), (1:3)))) # yields tuple of index ranges
(-2:5, 1:3)

julia> as(CartesianIndices, ((-2:5), (1:3))) # calls constructor
CartesianIndices((-2:5, 1:3))

julia> as(String, :hello) # converts symbol to string
"hello"

julia> as(Symbol, "hello") # converts string to symbol
:hello
```

Another usage is:

``` julia
as(T)
```

which yields a callable object that converts its argument to type `T`. This can be useful
with `map`. For instance:

``` julia
map(as(Int), dims)
```

to convert `dims` to a tuple (or array) of `Int`s.

Additional conversions becomes possible if another package such as
[`TwoDimensonal`](https://github.com/emmt/TwoDimensional.jl) is loaded.

For more general lazy conversions, consider [`lazymap` in this package](#lazy-maps).


## Lazy maps

Given an array or an iterator `A`, a *lazy map* `B` is built by one of:

```julia
B = lazymap([T::Type,] f, A)
B = lazymap([T::Type,] f, A::AbstractArray, f_inv = inverse(f))
```

which yield a view `B` of the array or iterator `A` such that the `i`-th element of `B` is
given by `Bᵢ = as(T, f(Aᵢ))` with `Aᵢ` the `i`-th element of `A` and using the [`as`
method](#cast-value-to-type) to convert the result of `f(Aᵢ)` to type `T`.

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

```julia
C = lazymap(T::Type, A)
```

builds an object `C` that lazily maps the **constructor** `T` to the elements of `A`. This
is not exactly the same as:

```julia
B = lazymap(T::Type, identity, A)
```

which builds an object `B` that lazily **converts** the elements of `A` to type `T`. In
other words, the `i`-th element of `C` is given by `Cᵢ = T(Aᵢ)::T`, while the `i`-th element
of `B` is given by `Bᵢ = as(T, Aᵢ)`. In both cases, it is asserted that `Cᵢ` and `Bᵢ` are of
type `T`. The two are equivalent if `T` is a numeric type (a sub-type of `Number`).


## Unitless value or type

`unitless(x)` yields `x` without its units, if any. `x` can be a number or a numeric type.
In the latter case, `unitless` behaves like `bare_type` described below.


## Bare, real, and floating-point numerical types

`bare_type(x)` yields the bare numeric type of `x` (a numeric value or type). If this
method is not extended for a specific type, the fallback implementation yields
`typeof(one(x))`. With more than one argument, `bare_type(args...)` yields the type
resulting from promoting the bare numeric types of `args...`. With no argument,
`bare_type()` yields `TypeUtils.BareNumber` the union of bare numeric types that may be
returned by this method.

`real_type(x)` yields the bare real type of `x` (a numeric value or type). If this method
is not extended for a specific type, the fallback implementation yields
`typeof(one(real(x))`. With more than one argument, `real_type(args...)` yields the type
resulting from promoting the bare real types of `args...`. With no argument, `real_type()`
yields `Real` the super-type of types that may be returned by this method.

The only difference between `bare_type` and `real_type` is how they treat complex numbers.
The former preserves the complex kind of its argument while the latter always returns a
real type. You may assume that `real_type(x) = real(bare_type(x))`. Conversely,
`convert_bare_type(T,x)` yields a complex result if `T` is complex and a real result if
`T` is real whatever `x`, while `convert_real_type(T,x)` yields a complex result if `x` is
complex and a real result if `x` is real, only the real part of `T` matters for
`convert_real_type(T,x)`. See examples below.

`floating_point_type(args...)` yields a floating-point type appropriate to represent the
bare real type of `args...`. With no argument, `floating_point_type()` yields
`AbstractFloat` the super-type of types that may be returned by this method. You may
consider `floating_point_type(args...)` as an equivalent to to`float(real_type(args...))`.
The floating-point type can be seen as the numerical precision for computations involving
`args...`.

`convert_bare_type(T,x)` converts the bare numeric type of `x` to the bare numeric type of
`T` while preserving the units of `x` if any. Argument `x` may be a number or a numeric
type, while argument `T` must be a numeric type. If `x` is one of `missing`, `nothing`,
`undef`, or the type of one of these singletons, `x` is returned.

`convert_real_type(T,x)` converts the bare real type of `x` to the bare real type of `T`
while preserving the units of `x` if any. Argument `x` may be a number or a numeric type,
while argument `T` must be a numeric type. If `x` is one of `missing`, `nothing`, `undef`,
or the type of one of these singletons, `x` is returned.

`convert_floating_point_type(T,x)` converts the bare real type of `x` to the suitable
floating-point type for type `T` while preserving the units of `x` if any. Argument `x`
may be a number or a numeric type, while argument `T` must be a numeric type. If `x` is
one of `missing`, `nothing`, `undef`, or the type of one of these singletons, `x` is
returned. You may consider `convert_floating_point_type(T,x)` as an equivalent to to
`convert_real_type(float(real_type(T)),x)`.

## Numerical precision

`get_precision(x)` or `get_precision(typeof(x))` yield the numerical precision of `x`. The
precision is a dimensionless floating-point type. If `x` is a floating-point value, the
precision of `x` is its floating-point type; if `x` stores floating-point values, the
precision is their promoted floating-point type; otherwise, the precision is
`AbstractFloat`.

To adapt the numerical precision of an object `x`, call:

``` julia
adapt_precision(T::Type{<:AbstractFloat}, x)
```

which yields an object similar to `x` but with numerical precision specified by the
floating-point type `T`. If `x` has already the required precision or if setting its
precision is irrelevant or not implemented, `x` is returned unchanged. Setting the
precision shall not change the dimensions of dimensionful numbers. If `T` is
`AbstractFloat`, the default floating-point type `Float64` is assumed.

For a number `x`, `adapt_precision(T, x)` behaves as `convert_real_type(T, x)` and
`adapt_precision(T, typeof(x))` may be used to infer the corresponding type with precision
`T`.

Example:

```julia
julia> adapt_precision(Float32, (1, 0x07, ("hello", 1.0, 3.0 - 2.0im, π)))
(1.0f0, 7.0f0, ("hello", 1.0f0, 3.0f0 - 2.0f0im, 3.1415927f0))
```

Basically, `adapt_precision` supports numbers, arrays of numbers, matrix factorizations,
and (named) tuples. It can be specialized for other object types defined in foreign
packages by specializing:

```julia
TypeUtils.adapt_precision(::Type{T}, x::SomeType) where {T<:TypeUtils.Precision} = ...
```

where `SomeType` is the object type and `TypeUtils.Precision` is a public but non-exported
symbol defined in `TypeUtils` as the union of the concrete floating-point types of Julia:

```julia
const Precision = Union{Float16,Float32,Float64,BigFloat}
```

The restriction `T<:Precision` in the function signature is to make sure that this version
of `adapt_precision` is only called with a concrete floating-point type `T`.


## Parameter-less type

The call:

``` julia
parameterless(T)
```

yields the type `T` without parameter specifications. For example:

```julia
julia> parameterless(Vector{Float32})
Array
```

## Deal with array shape, size, and axes

The `TypeUtils` package provides the following types for array shape, size, and axes:

- `ArrayAxis` is an alias to the possible canonical types representing a single array
  index range, that is `AbstractUnitRange{eltype(Dims)}`.

- `ArrayAxes{N}` is an alias to the possible canonical types representing `N`-dimensional
  array axes, that is `NTuple{N,ArrayAxis}`. For any `N`-dimensional array `A`, `axes(A)
  <: ArrayAxes{N}` should hold.

- `ArrayShape{N}` is an alias to the `N`-tuple of array dimensions and/or index ranges to
  which `as_array_shape`, `as_array_size`, or `as_array_axes` are applicable (see below).

Note that `Dims{N}` is the same as `NTuple{N,Int}` in Julia and represents the result of
`size(A)` for the `N`-dimensional array `A`, so `eltype(Dims) = Int`.

Given a list `inds...` of array dimension lengths (integers) and/or index ranges
(integer-valued unit ranges), the following methods from `TypeUtils` are applicable:

- `as_array_shape(inds...)` yields canonical array axes (if `inds...` contains any index
  range other than `Base.OneTo`) or canonical array size (otherwise). The former is an
  instance of `ArrayAxes{N}`, the latter is an instance of `Dims{N}`.

- `as_array_axes(inds...)` yields canonical array axes, that is a `N`-tuple of type
  `ArrayAxes{N}`.

- `as_array_size(inds...)` yields canonical array size, that is a `N`-tuple of type
  `Dims{N}`.

Of course, these methods also accept that their arguments be specified by a tuple of
array dimension lengths and/or index ranges, that is `(inds...,)` instead of `inds...` in
the above example.

To deal with individual array dimension length and/or index range:

- `as_array_axis` converts its argument to a single array axis of type `ArrayAxis`.

- `as_array_dim` converts its argument to a single array dimension length of type
  `eltype(Dims)`, that is an `Int`.

The method `new_array(T, inds...)` yields a new array with elements of type `T` and shape
defined by `inds...`. Following the semantic of `as_array_shape`, the returned array is an
`OffsetArray{T}` if `inds...` contains any index range other than `Base.OneTo` and an
`Array{T}` otherwise. In the former case, an exception is thrown if the package
`OffsetArrays` has not been loaded.


## Deal with array element types

The `TypeUtils` package provides a few methods to deal with array element types:

* `promote_eltype(args...)` yields the promoted element type of the arguments `args...`
  which may be anything implementing the `eltype` method.

* `convert_eltype(T,A)` yields an object similar to `A` except that its elements have type
  `T`.

* `as_eltype(T,A)` yields an array which lazily converts its entries to type `T`. This can
  be seen as a memory-less version of `convert_eltype(T,A)`. The method `as_eltype` is
  similar to the method `of_eltype` provided by the
  [`MappedArrays`](https://github.com/JuliaArrays/MappedArrays.jl/tree/master) package.

Methods `convert_eltype(T,A)` and `as_eltype(T,A)` just return `A` itself if its elements
are of type `T`.


## Type of result returned by a function

The call:

``` julia
g = as_return(T, f)
```

yields a callable object such that `g(args...; kwds...)` lazily converts the value
returned by `f(args...; kwds...)` to the type `T`. Methods `return_type(g)` and
`parent(g)` can be used to respectively retrieve the type `T` and the original function
`f`. A similar kind of object be built with the composition operator:

``` julia
g = as(T)∘f
```

The method `return_type` may also be used as:

``` julia
T = return_type(f, argtypes...)
```

to infer the type `T` of the result returned by `f` when called with arguments of types
`argtypes...`.


## Destructure an object as values and restructure values as an object

It is sometime useful to collect the values stored by a structured object into simple
collection of values (a tuple or a vector). The reverse operation is also needed. The
call:

``` julia
vals = destructure(obj)
```

yields a tuple, `vals`, of the values of the structured object `obj` and, conversely, the
call:

``` julia
obj = restructure(T, vals)
```

builds a structured object of type `T` given `vals`, a tuple or a vector of its values.

For an immutable concrete object `obj`, the following identity holds:

``` julia
restructure(typeof(obj), destructure(obj)) === obj
```

It is also possible to destructure an object into a given vector of values:

``` julia
destructure!(vals, obj)
```

Optionally, in `restructure` and `destructure!` methods, keyword `offset` may be specified
to not start with the first value in `vals`.

Method `struct_length` yields the minimal number of values needed to destructure an
object.


## Compatibility with `Unitful`

The following examples illustrate the result of the methods provided by `TypeUtils`, first
with bare numbers and bare numeric types, then with quantities:

```julia
julia> using TypeUtils

julia> map(unitless, (2.1, Float64, true, ComplexF32))
(2.1, Float64, true, ComplexF32)

julia> map(bare_type, (1, 3.14f0, true, 1//3, π, 1.0 - 2.0im))
(Int64, Float32, Bool, Rational{Int64}, Irrational{:π}, Complex{Float64})

julia> map(real_type, (1, 3.14f0, true, 1//3, π, 1.0 - 2.0im))
(Int64, Float32, Bool, Rational{Int64}, Irrational{:π}, Float64)

julia> map(x -> convert_bare_type(Float32, x), (2, 1 - 0im, 1//2, Bool, Complex{Float64}))
(2.0f0, 1.0f0, 0.5f0, Float32, Float32)

julia> map(x -> convert_real_type(Float32, x), (2, 1 - 0im, 1//2, Bool, Complex{Float64}))
(2.0f0, 1.0f0 + 0.0f0im, 0.5f0, Float32, ComplexF32)

julia> using Unitful

julia> map(unitless, (u"2.1GHz", typeof(u"2.1GHz")))
(2.1, Float64)

julia> map(bare_type, (u"3.2km/s", u"5GHz", typeof((0+1im)*u"Hz")))
(Float64, Int64, Complex{Int64})

julia> map(real_type, (u"3.2km/s", u"5GHz", typeof((0+1im)*u"Hz")))
(Float64, Int64, Int64)
```


## Rationale

The following example shows a first attempt to use `bare_type` to implement efficient
in-place multiplication of an array (whose element may have units) by a real factor (which
must be dimensionless in this context):

```julia
function scale!(A::AbstractArray, α::Number)
    alpha = convert_bare_type(eltype(A), α)
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end
```

An improvement is to realize that when `α` is a real while the entries of `A` are
complexes, it is more efficient to multiply the entries of `A` by a real-valued multiplier
rather than by a complex one. Implementing this is as simple as replacing
`convert_bare_type` by `convert_real_type` to only convert the bare real type of the
multiplier while preserving its complex/real kind:

```julia
function scale!(A::AbstractArray, α::Number)
    alpha = convert_real_type(eltype(A), α)
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end
```

This latter version consistently and efficiently deals with `α` being real while the
entries of `A` are reals or complexes, and with `α` and the entries of `A` being
complexes. If `α` is a complex and the entries of `A` are reals, the statement `A[i] *=
alpha` will throw an `InexactConversion` if the imaginary part of `α` is not zero. This
check is probably optimized out of the loop by Julia but, to handle this with guaranteed
no loss of efficiency, the code can be written as:

```julia
function scale!(A::AbstractArray, α::Union{Real,Complex})
    alpha = if α isa Complex && bare_type(eltype(A)) isa Real
        convert(real_type(eltype(A)), α)
    else
        convert_real_type(eltype(A), α)
    end
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end
```

The restriction `α::Union{Real,Complex}` accounts for the fact that in-place
multiplication imposes a dimensionless multiplier. Since the test leading to the
expression used for `alpha` is based on the types of the arguments, the branch is
eliminated at compile time and the type of `alpha` is known by the compiler. The
`InexactConversion` exception may then only be thrown by the call to `convert` in the
first branch of the test.

This seemingly very specific case was in fact the key point to allow for packages such as
[LazyAlgebra](https://github.com/emmt/LazyAlgebra.jl) or
[LinearInterpolators](https://github.com/emmt/LinearInterpolators.jl) to work seamlessly
on arrays whose entries may have units. The `TypeUtils` (formerly `Unitless`) package was
created to cover this need as transparently as possible.


## Related things

There exist objects or packages with functionalities similar to those provided by
`TypeUtils` but with differences that justify the existence of this package.

* Compared to `Iterators.map(f, A)` which is always an iterator, the object returned by
  `lazymap(f, A)` is an (abstract) array if `A` is an array, an iterator otherwise.

* `mapview(f, A)` using [`FlexiMaps`](https://github.com/JuliaAPlavin/FlexiMaps.jl) is very
  similar to `lazymap(f, A)` but, with `FlexiMaps`, there is no equivalent to `lazymap(T, f,
  A)` to specify an element type for the mapped view and directly build a
  `FlexiMaps.MappedArray` with a given element type yields an abstract array whose `eltype`
  is not the type of the result returned by `getindex` while `lazymap` takes care of
  converting this result correctly. In `FlexiMaps`, there is no shortcut to build a
  read-only view if the inverse function is known.

* `mappedarray(f, A)` and `mappedarray(f, inv_f, A)` using
  [`MappedArrays`](https://github.com/JuliaArrays/MappedArrays.jl) are similar to
  `lazymap(f, A)` and `lazymap(f, A, inv_f)` but:
  - `f_inv` is not automatically inferred by `mappedarray` if not specified;
  - `A` must be an array in `mappedarray`;
  - the element type of a mapped array is inferred and cannot be specified as for a lazy map
    which is type-stable with respect to its `eltype`.

* `lazymap(T, A)` is the analogue of `of_eltype(T, A)` and `as_eltype(T, A)` respectively
  using [`MappedArrays`](https://github.com/JuliaArrays/MappedArrays.jl) and
  [`TypeUtils`](https://github.com/emmt/TypeUtils.jl).

* `LazyMaps` does not implement lazily mapping multiple arrays, a possibility offered by
  `MappedArrays`, but this may be emulated by combining `LazyMaps` and
  [`ZippedArrays`](https://github.com/emmt/ZippedArrays.jl).

* `BroadcastArray(f, A)` and `BroadcastArray{T}(f, A)` using
  [`LazyArrays`](https://github.com/JuliaArrays/LazyArrays.jl) is similar to `lazymap(f, A)`
  and `lazymap(T, f, A)` but broadcast arrays are read-only and restricted to arrays while
  lazy maps are writable if inverse function is known or specified and can be used over
  other collections than arrays.

As shown by [benchmark tests](./test/benchmarks.jl) for `LazyMaps` and these different
packages, evaluating `B[i]` for an object `B` lazily representing `f.(A)` is as fast as
calling `f(A[i])`. Also any of these objects can be used with no allocations and, except
`BroadcastArray`, no construction overheads compared to `f(A[i])`. A `BroadcastArray` using
[`LazyArrays`](https://github.com/JuliaArrays/LazyArrays.jl) is as fast as `mapreduce` for
reductions (like a `sum`) of the broadcast array which provides some speedup for large
arrays.

Direct dependencies:

* [`InverseFunctions`](https://github.com/JuliaMath/InverseFunctions.jl) is used to infer
  inverse functions.


## Installation

`TypeUtils` can be installed as any other official Julia packages. For example:

```julia
using Pkg
Pkg.add("TypeUtils")
```
