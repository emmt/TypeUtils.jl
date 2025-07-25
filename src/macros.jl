"""
    TypeUtils.@public args...

declares `args...` as being `public` even though they are not exported. For Julia version <
1.11, this macro does nothing. Using this macro also avoid errors with CI and coverage
tools.

"""
macro public(args::Union{Symbol,Expr}...)
    VERSION ≥ v"1.11.0-DEV.469" ? esc(Expr(:public, map(
        x -> x isa Symbol ? x :
            x isa Expr && x.head == :macrocall ? x.args[1] :
            error("unexpected argument `$x` to `@public`"), args)...)) : nothing
end
VERSION ≥ v"1.11.0-DEV.469" && @public @public

"""
    @assert_floating_point A B ...

throws an `ArgumentError` exception if any of the variables `A`, `B`, etc. does not use
floating-point to store its value(s).

See also [`assert_floating_point`](@ref).

"""
macro assert_floating_point(args::Symbol...)
    code = [:(assert_floating_point($(QuoteNode(arg)), $(esc(arg)))) for arg in args]
    return quote
        $(code...)
    end
end
