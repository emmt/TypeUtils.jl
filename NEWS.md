# User visible changes in `TypeUtils`

# Version 1.2.0

- `to_same_type(x...)` is a substitute to `promote(x...)` that warrants that returned
  instances have the same type and that calls `as(T,x)`, not `convert(T,x)`, if any
  conversion is needed.

- `to_same_concrete_type(T...)` is a substitute to `promote_type(T...)` that throws an
  exception if types `T...` cannot be promoted to a common concrete type.

# Version 1.1.0

- `convert_eltype` can be applied to a number.

# Version 1.0.0

- All methods and types formerly provided by
  [`Unitless`](https://github.com/emmt/Unitless.jl) are now provided by
  `TypeUtils` which supersedes `Unitless`.

# Version 0.3.8

- `TwoDimensional` is no longer an extension because `TwoDimensional` version
  0.5 directly extends `convert` as expected by `TypeUtils`.

# Version 0.3.7

- Method `as_return(T, f)` builds an instance of `TypeStableFunction{T}` which
  is a sub-type of `AbstractTypeStableFunction{T}`. These two types are both
  exported.

# Version 0.3.6

- Methods `destructure`, `destructure!`, and `restructure` are inline.

# Version 0.3.5

- Methods `destructure`, `destructure!`, `restructure`, and `struct_length`
  can deal with tuples.

# Version 0.3.4

- New methods `vals = destructure(obj)` or `destructure!(vals, obj)`, and `obj
  = restructure(T, vals)` to destructure and object `obj` as a tuple or vector
  of its values and, conversely, to rebuild an object of type `T` from its
  values.

- New method `struct_length` yields the number of values needed to destructure
  an object.

# Version 0.3.3

- Fix `convert_eltype(T,A)` when `A` is a range.

- `parameterless` is now implemented as `constructorof` in
  [`ConstructionBase`](https://github.com/JuliaObjects/ConstructionBase.jl).


# Version 0.3.2

- `convert_eltype(T,A)` yields a tuple if `A` is a tuple.


# Version 0.3.1

- `convert_eltype(T,A)` yields a range if `A` is a range.


# Version 0.3.0

This is the first version as an official Julia package.


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
