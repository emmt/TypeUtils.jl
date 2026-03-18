module TypeUtilsLinearAlgebraExt

if isdefined(Base, :get_extension)
    using TypeUtils, LinearAlgebra
else
    using ..TypeUtils, ..LinearAlgebra
end

# Convert element type for LinearAlgebra factorizations.
# `LinearAlgebra.Factorization{T}(A)` can be used to convert element type of `A` for QR,
# LinearAlgebra.QRCompactWY, QRPivoted, LQ, Cholesky, CholeskyPivoted, LU, LDLt,
# BunchKaufman, SVD, etc.
TypeUtils.convert_eltype(::Type{T}, A::Factorization{T}) where {T} = A
TypeUtils.convert_eltype(::Type{T}, A::Factorization) where {T} = Factorization{T}(A)
if VERSION < v"1.7.0-beta1"
    # For old Julia versions, the above is not sufficient for SVD.
    TypeUtils.convert_eltype(::Type{T}, A::SVD{T}) where {T} = A
    TypeUtils.convert_eltype(::Type{T}, A::SVD) where {T} =
        SVD(TypeUtils.convert_eltype(T, A.U), TypeUtils.convert_eltype(real(T), A.S), TypeUtils.convert_eltype(T, A.Vt))
end
TypeUtils.convert_eltype(::Type{T}, A::Hessenberg{T}) where {T} = A
TypeUtils.convert_eltype(::Type{T}, A::Hessenberg) where {T} = throw(
    ArgumentError(
        "changing element type of Hessenberg decomposition is not supported, consider re-computing the decomposition"
    )
)

# For `Adjoint` and `Transpose`, we want to preserve this structure.
for S in (:Adjoint, :Transpose)
    @eval begin
        TypeUtils.convert_eltype(::Type{T}, A::$S{T}) where {T} = A
        TypeUtils.convert_eltype(::Type{T}, A::$S) where {T} = $S(TypeUtils.convert_eltype(T, parent(A)))
    end
end

# Get and adapt precision for LinearAlgebra factorizations.
TypeUtils.get_precision(::Type{<:Factorization{T}}) where {T} = TypeUtils.get_precision(T)
TypeUtils.adapt_precision(::Type{T}, A::Factorization{T}) where {T <: TypeUtils.Precision} = A
TypeUtils.adapt_precision(::Type{T}, A::Factorization{S}) where {T <: TypeUtils.Precision, S} =
    TypeUtils.convert_eltype(TypeUtils.adapt_precision(T, S), A)

end # module
