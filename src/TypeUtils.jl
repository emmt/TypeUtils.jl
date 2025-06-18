module TypeUtils

export
    @assert_floating_point,
    ArrayAxes,
    ArrayAxis,
    ArrayShape,
    RelaxedArrayShape,
    AbstractTypeStableFunction,
    TypeStableFunction,
    as,
    as_array_axes,
    as_array_axis,
    as_array_dim,
    as_array_shape,
    as_array_size,
    as_eltype,
    as_return,
    assert_floating_point,
    bare_type,
    convert_bare_type,
    convert_eltype,
    convert_floating_point_type,
    convert_real_type,
    destructure!,
    destructure,
    floating_point_type,
    is_signed,
    nearest,
    new_array,
    parameterless,
    promote_eltype,
    real_type,
    restructure,
    return_type,
    struct_length,
    to_same_type,
    to_same_concrete_type,
    unitless

using Base: OneTo

if !isdefined(Base, :get_extension)
    using Requires
end

include("macros.jl")

@public BareNumber
@public Converter
@public Dim
@public Unsupported

include("types.jl")
include("methods.jl")
include("numbers.jl")
include("arrays.jl")
include("funcs.jl")
include("structs.jl")

function __init__()
    @static if !isdefined(Base, :get_extension)
        # Extend methods when other packages are loaded.
        @require Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d" include(
            "../ext/TypeUtilsUnitfulExt.jl")
        @require OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881" include(
            "../ext/TypeUtilsOffsetArraysExt.jl")
    end
end

end # module TypeUtils
