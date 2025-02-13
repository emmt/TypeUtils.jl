# User visible changes in `TypeUtils`

- New `nearest(T,x)` method to return the value of type `T` that is the nearest to `x`.
  For `T` integer and `x` real, it can be seen as rounding with clamping to avoid
  overflows.
- Non-exported `TypeUtils.@public` macro is now public.

# Version 1.4.1

- `AbstractUnitRange{<:Integer}` has been replaced by `AbstractRange{<:Integer}` as the
  eligible type for specifying an array axis. Methods `as_array_axis`, `as_array_axes`,
  and `as_array_shape` convert index ranges to `AbstractUnitRange{Int}` throwing an
  exception if any range does not have unit step. The downside is that this must be
  explicitly checked for non-`AbstractUnitRange{<:Integer}` ranges while it can be
  inferred by type (it is a trait) for `AbstractUnitRange{<:Integer}` ranges.

- The exported alias `RelaxedArrayShape{N}` reflects the above change. Exported alias
  `ArrayShape{N}` is restricted to integers and integer-valued unit-ranges which better
  correspond to Julia's way of representing an `N`-dimensional array shape.

- New macro `TypeUtils.@public` used to declare non-exported public symbols. Does nothing
  for Julia versions older than 1.11.

- For other packages, it may be sufficient to extend `convert_eltype(T, X)` to a given
  type `X` to have `convert_eltype(T, x::X)` work for instances `x` of that type.

- `TypeStableFunction(f, argtypes...)` better tries to infer a suitable concrete type for
  `f` with arguments of types `argtypes...`.

# Version 1.4.0

- Non-exported type `TypeUtils.Unsupported` may be used to provide a fallback
  implementation of a method for given types of arguments that is only supported when some
  extension is loaded.

- `new_array(T, inds...)` creates an array with elements of type `T` and shape defined by
  `inds...`. The returned array is an `OffsetArray{T}` if `inds...` contains any index
  range other than `Base.OneTo` and an `Array{T}` otherwise. In the former case, an
  exception is thrown if the package `OffsetArrays` has not been loaded.

# Version 1.3.0

Add a few types and methods related to array size and axes:

- `ArrayAxis = AbstractUnitRange{eltype(Dims)}` is an alias to the possible canonical
  types representing a single array index range.

- `ArrayAxes{N} = NTuple{N,ArrayAxis}` is an alias to the possible canonical types
  representing `N`-dimensional array axes.

- `ArrayShape{N}` is an alias to the `N`-tuple of array dimensions and/or index ranges
  to which `as_array_shape`, `as_array_size`, or `as_array_axes` are applicable.

- `as_array_shape` converts its argument(s) to canonical array axes or to canonical array
  size.

- `as_array_axes` converts its argument(s) to canonical array axes, that is a `N`-tuple of
  type `ArrayAxes{N}`. `as_array_axis` converts its argument to a single array axis
  of type `ArrayAxis`.

- `as_array_size` converts its argument(s) to a canonical array size, that is a `N`-tuple
  of type `Dims{N}`. `as_array_dim` converts its argument to a single array dimension
  length of type `eltype(Dims)`.

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
