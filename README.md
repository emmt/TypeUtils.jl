# Casting to type in Julia

[![Build Status](https://github.com/emmt/AsType.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/AsType.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/AsType.jl?svg=true)](https://ci.appveyor.com/project/emmt/AsType-jl) [![Coverage](https://codecov.io/gh/emmt/AsType.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/AsType.jl)

`AsType` is a small [Julia](https://www.julialang.org) package providing a
single method, `as`, designed to *cast* an argument to a given type. The name
was inspired by the built-in Zig function
[`@as`](https://ziglang.org/documentation/master/#as).

A first usage is:

``` julia
as(T, x)
```

which yields `x` converted to type `T`. This behaves like a lazy version of
`convert(T,x)::T` doing nothing if `x` is already of type `T` and performing
the conversion and the type assertion otherwise.

The `as` method calls `convert` by default but implements a number of
conversions not supported by `convert`. The `as` method is therefore a bit more
versatile than `convert` while relaxing the bother to remember which function
or constructor to call to perform the intended conversion. For example:

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
