# User visible changes in `TypeUtils`

# Version 0.2.4

- Fix a typo causing `as_return` to fail when applied to an `AsReturn` object.


# Version 0.2.3

- New method `as_return(T, f)` yields a callable object that behaves like `f`
  when called as a function except that it lazily converts the value returned
  by `f` to type `T`.

- New method `return_type(f, argtypes...)` to infer the type of the result
  returned by the callable object `f` when called with arguments of types
  `argtypes...`.


# Version 0.2.2

- New method `promote_eltype(args...)` to yield the promoted element type of
  `args...`.

- New method `convert_eltype(T,A)` to convert the element type of `A` to be `T`.

- New method `as_eltype(T,A)` to lazily convert the element type of `A` to be `T`.


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
