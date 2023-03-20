# Tiny utility for casting to type in Julia

[![Build Status](https://github.com/emmt/AsType.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/emmt/AsType.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Build Status](https://ci.appveyor.com/api/projects/status/github/emmt/AsType.jl?svg=true)](https://ci.appveyor.com/project/emmt/AsType-jl) [![Coverage](https://codecov.io/gh/emmt/AsType.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/emmt/AsType.jl)

`AsType` is a super tiny [Julia](https://www.julialang.org) package providing a
single method, `as`, designed to *cast* its argument to a given type. The name
was inspired by the Zig built-in function
[`@as`](https://ziglang.org/documentation/master/#as).

A first possible usage is:

``` julia
as(T, x)
```

yields `x` converted to type `T`. It behaves like a lazy version of
`convert(T,x)::T` doing nothing if `x` is already of type `T` and performing
the conversion and the type assertion otherwise.

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
