# User visible changes in `TypeUtils`

# Version 0.2.1

- Method `parameterless(T)` to get the type `T` without parameter
  specifications. For example:

  ```julia
  julia> parameterless(Vector{Float32})
  Array
  ```


# Version 0.2.0

- Methods `as(T,x)` to convert `x` to type `T`.

- Call `f = as(T)` to build a callable object such that `f(x)` is the same as
  `as(T,x)`.

- Extension for [`TwoDimensonal`](https://github.com/emmt/TwoDimensional.jl).

- Package was previously named `AsType` and `CastType`.
