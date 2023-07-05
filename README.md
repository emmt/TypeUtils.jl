# Dealing with types in Julia

[![Build Status](https://github.com/emmt/TypeUtils.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/TypeUtils.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/TypeUtils.jl?svg=true)](https://ci.appveyor.com/project/emmt/TypeUtils-jl) [![Coverage](https://codecov.io/gh/emmt/TypeUtils.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/TypeUtils.jl)

Package `TypeUtils` provides useful methods to deal with types in
[Julia](https://www.julialang.org).

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


## Parameter-less type

The method:

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

* `convert_eltype(T,A)` yields an array identical to `A` except that its
  elements have type `T`.

* `as_eltype(T,A)` yields an array which lazily converts its entries to type
  `T`. This can be seen as a memory-less version of `convert_eltype(T,A)`. The
  method `as_eltype` is similar to the method `of_eltype` provided by the
  [`MappedArrays`](https://github.com/JuliaArrays/MappedArrays.jl/tree/master)
  package.

Methods `convert_eltype(T,A)` and `as_eltype(T,A)` just return `A` itself if
its elements are of type `T`.
