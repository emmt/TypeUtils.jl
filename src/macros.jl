"""
    TypeUtils.@public a b c ...
    TypeUtils.@public a, b, c, ...

Declare symbols `a`, `b`, `c`, etc. as being `public` even though they are not exported. For
Julia version < 1.11, this macro does nothing. Using this macro also avoid errors with CI
and coverage tools.

"""
macro public(args...)
    VERSION ≥ v"1.11.0-DEV.469" || return nothing
    # `@public a b c` and `@public(a, b, c)` are the same, but `@public a, b, c` is
    # different. Make `xs` a tuple or a vector of public symbols in these different cases.
    xs = args isa Tuple{Expr} && args[1].head === :tuple ? args[1].args : args
    return esc(Expr(:public, map(
        x -> x isa Symbol ? x : x isa Expr && x.head == :macrocall ? x.args[1] :
            error("unexpected argument `$x` to `@public`"), xs)...))
end

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
