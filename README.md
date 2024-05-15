# Dealing with types in Julia

[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](./LICENSE.md)  [![Build Status](https://github.com/emmt/TypeUtils.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/TypeUtils.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/TypeUtils.jl?svg=true)](https://ci.appveyor.com/project/emmt/TypeUtils-jl) [![Coverage](https://codecov.io/gh/emmt/TypeUtils.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/TypeUtils.jl)

Package `TypeUtils` provides useful methods to deal with types in
[Julia](https://www.julialang.org) and facilitate coding with numbers whether
they have units or not. The package provides methods to strip units from
numbers or numeric types, convert the numeric type of quantities (not their
units), determine appropriate numeric type to carry computations mixing numbers
with different types and/or units. These methods make it easy to write code
that works consistently for numbers with any units (including none). The
intention is that the `TypeUtils` package automatically extends its exported
methods when packages such as
[`Unitful`](https://github.com/PainterQubits/Unitful.jl) are loaded.


## Cast value to type

The method, `as` is designed to *cast* a value to a given type. The name was
inspired by the built-in [Zig](https://ziglang.org/) function
[`@as`](https://ziglang.org/documentation/master/#as).

A first usage is:

``` julia
as(T, x)
```

which yields `x` converted to type `T`. This behaves like a lazy version of
`convert(T,x)::T` doing nothing if `x` is already of type `T` and performing
the conversion and the type assertion otherwise.

By default, the `as` method calls `convert` only if needed but also implements
a number of conversions not supported by `convert`. The `as` method is
therefore a bit more versatile than `convert` while relaxing the bother to
remember which function or constructor to call to efficiently perform the
intended conversion. For example:

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

which yields a callable object that converts its argument to type `T`. This can
be useful with `map`. For instance:

``` julia
map(as(Int), dims)
```

to convert `dims` to a tuple (or array) of `Int`s.

Additional conversions becomes possible if another package such as
[`TwoDimensonal`](https://github.com/emmt/TwoDimensional.jl) is loaded.


## Unitless value or type

`unitless(x)` yields `x` without its units, if any. `x` can be a number or a
numeric type. In the latter case, `unitless` behaves like `bare_type` described
below.


## Bare, real, and floating-point numerical types

`bare_type(x)` yields the bare numeric type of `x` (a numeric value or type).
If this method is not extended for a specific type, the fallback implementation
yields `typeof(one(x))`. With more than one argument, `bare_type(args...)`
yields the type resulting from promoting the bare numeric types of `args...`.
With no argument, `bare_type()` yields `TypeUtils.BareNumber` the union of bare
numeric types that may be returned by this method.

`real_type(x)` yields the bare real type of `x` (a numeric value or type). If
this method is not extended for a specific type, the fallback implementation
yields `typeof(one(real(x))`. With more than one argument, `real_type(args...)`
yields the type resulting from promoting the bare real types of `args...`. With
no argument, `real_type()` yields `Real` the super-type of types that may be
returned by this method.

The only difference between `bare_type` and `real_type` is how they treat
complex numbers. The former preserves the complex kind of its argument while
the latter always returns a real type. You may assume that `real_type(x) =
real(bare_type(x))`. Conversely, `convert_bare_type(T,x)` yields a complex
result if `T` is complex and a real result if `T` is real whatever `x`, while
`convert_real_type(T,x)` yields a complex result if `x` is complex and a real
result if `x` is real, only the real part of `T` matters for
`convert_real_type(T,x)`. See examples below.

`floating_point_type(args...)` yields a floating-point type appropriate to
represent the bare real type of `args...`. With no argument,
`floating_point_type()` yields `AbstractFloat` the super-type of types that may
be returned by this method. You may consider `floating_point_type(args...)` as
an equivalent to to`float(real_type(args...))`. The floating-point type can be
seen as the numerical precision for computations involving `args...`.

`convert_bare_type(T,x)` converts the bare numeric type of `x` to the bare
numeric type of `T` while preserving the units of `x` if any. Argument `x`
may be a number or a numeric type, while argument `T` must be a numeric type.
If `x` is one of `missing`, `nothing`, `undef`, or the type of one of these
singletons, `x` is returned.

`convert_real_type(T,x)` converts the bare real type of `x` to the bare real
type of `T` while preserving the units of `x` if any. Argument `x` may be a
number or a numeric type, while argument `T` must be a numeric type. If `x`
is one of `missing`, `nothing`, `undef`, or the type of one of these
singletons, `x` is returned.

`convert_floating_point_type(T,x)` converts the bare real type of `x` to the
suitable floating-point type for type `T` while preserving the units of `x`
if any. Argument `x` may be a number or a numeric type, while argument `T`
must be a numeric type. If `x` is one of `missing`, `nothing`, `undef`, or
the type of one of these singletons, `x` is returned. You may consider
`convert_floating_point_type(T,x)` as an equivalent to
to `convert_real_type(float(real_type(T)),x)`.

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


## Deal with array element types

The `TypeUtils` package provides a few methods to deal with array element
types:

* `promote_eltype(args...)` yields the promoted element type of the arguments
  `args...` which may be anything implementing the `eltype` method.

* `convert_eltype(T,A)` yields an object similar to `A` except that its
  elements have type `T`.

* `as_eltype(T,A)` yields an array which lazily converts its entries to type
  `T`. This can be seen as a memory-less version of `convert_eltype(T,A)`. The
  method `as_eltype` is similar to the method `of_eltype` provided by the
  [`MappedArrays`](https://github.com/JuliaArrays/MappedArrays.jl/tree/master)
  package.

Methods `convert_eltype(T,A)` and `as_eltype(T,A)` just return `A` itself if
its elements are of type `T`.


## Type of result returned by a function

The call:

``` julia
g = as_return(T, f)
```

yields a callable object such that `g(args...; kwds...)` lazily converts the
value returned by `f(args...; kwds...)` to the type `T`. Methods
`return_type(g)` and `parent(g)` can be used to respectively retrieve the type
`T` and the original function `f`. A similar kind of object be built with the
composition operator:

``` julia
g = as(T)∘f
```

The method `return_type` may also be used as:

``` julia
T = return_type(f, argtypes...)
```

to infer the type `T` of the result returned by `f` when called with arguments
of types `argtypes...`.


## Destructure an object as values and restructure values as an object

It is sometime useful to collect the values stored by a structured object into
simple collection of values (a tuple or a vector). The reverse operation is
also needed. The call:

``` julia
vals = destructure(obj)
```

yields a tuple, `vals`, of the values of the structured object `obj` and,
conversely, the call:

``` julia
obj = restructure(T, vals)
```

builds a structured object of type `T` given `vals`, a tuple or a vector of its
values.

For an immutable concrete object `obj`, the following identity holds:

``` julia
restructure(typeof(obj), destructure(obj)) === obj
```

It is also possible to destructure an object into a given vector of values:

``` julia
destructure!(vals, obj)
```

Optionally, in `restructure` and `destructure!` methods, keyword `offset` may
be specified to not start with the first value in `vals`.

Method `struct_length` yields the minimal number of values needed to
destructure an object.


## Compatibility with `Unitful`

The following examples illustrate the result of the methods provided by
`TypeUtils`, first with bare numbers and bare numeric types, then with
quantities:

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

The following example shows a first attempt to use `bare_type` to implement
efficient in-place multiplication of an array (whose element may have units) by
a real factor (which must be unitless in this context):

```julia
function scale!(A::AbstractArray, α::Number)
    alpha = convert_bare_type(eltype(A), α)
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end
```

An improvement is to realize that when `α` is a real while the entries of `A`
are complexes, it is more efficient to multiply the entries of `A` by a
real-valued multiplier rather than by a complex one. Implementing this is as
simple as replacing `convert_bare_type` by `convert_real_type` to only convert
the bare real type of the multiplier while preserving its complex/real kind:

```julia
function scale!(A::AbstractArray, α::Number)
    alpha = convert_real_type(eltype(A), α)
    @inbounds @simd for i in eachindex(A)
        A[i] *= alpha
    end
    return A
end
```

This latter version consistently and efficiently deals with `α` being real
while the entries of `A` are reals or complexes, and with `α` and the entries
of `A` being complexes. If `α` is a complex and the entries of `A` are reals,
the statement `A[i] *= alpha` will throw an `InexactConversion` if the
imaginary part of `α` is not zero. This check is probably optimized out of the
loop by Julia but, to handle this with guaranteed no loss of efficiency, the
code can be written as:

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
multiplication imposes a unitless multiplier. Since the test leading to the
expression used for `alpha` is based on the types of the arguments, the branch
is eliminated at compile time and the type of `alpha` is known by the compiler.
The `InexactConversion` exception may then only be thrown by the call to
`convert` in the first branch of the test.

This seemingly very specific case was in fact the key point to allow for
packages such as [LazyAlgebra](https://github.com/emmt/LazyAlgebra.jl) or
[LinearInterpolators](https://github.com/emmt/LinearInterpolators.jl) to work
seamlessly on arrays whose entries may have units. The `TypeUtils` (formerly
`Unitless`) package was created to cover this need as transparently as
possible.


## Installation

`TypeUtils` can be installed as any other official Julia packages. For example:

```julia
using Pkg
Pkg.add("TypeUtils")
```
