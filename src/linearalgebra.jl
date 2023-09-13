using LinearAlgebra
using LinearAlgebra: MulAddMul, mul!, AdjOrTrans

@inline LinearAlgebra.generic_matvecmul!(C::OffsetVector, fA::Function, A::AbstractVecOrMat, B::AbstractVector,
                            alpha, beta) = unwrap_matvecmul!(C, fA, A, B, alpha, beta)

@inline function unwrap_matvecmul!(C::OffsetVector, fA, A::AbstractVecOrMat, B::AbstractVector,
                            alpha, beta)

    mB_axis = Base.axes1(B)
    mA_axis, nA_axis = axes(fA(A))

    if mB_axis != nA_axis
        throw(DimensionMismatch("mul! can't contract axis $(UnitRange(nA_axis)) from A with axes(B) == ($(UnitRange(mB_axis)),)"))
    end
    if mA_axis != Base.axes1(C)
        throw(DimensionMismatch("mul! got axes(C) == ($(UnitRange(Base.axes1(C))),), expected $(UnitRange(mA_axis))"))
    end

    mul!(no_offset_view(C), fA(no_offset_view(A)), no_offset_view(B), alpha, beta)
    C
end

# The signatures of these differs from LinearAlgebra's *only* on C:
@inline LinearAlgebra.generic_matmatmul!(C::OffsetMatrix, fA::Function, fB::Function, A::AbstractMatrix, B::AbstractMatrix,
                             alpha, beta) = unwrap_matmatmul!(C, fA, fB, A, B, alpha, beta)

@inline LinearAlgebra.generic_matmatmul!(C::Union{OffsetMatrix, OffsetVector}, fA::Function, fB::Function, A::AbstractVecOrMat, B::AbstractVecOrMat,
                             alpha, beta) = unwrap_matmatmul!(C, fA, fB, A, B, alpha, beta)

@inline LinearAlgebra.generic_matmatmul!(C::AdjOrTrans{<:Any, <:OffsetArray}, fA::Function, fB::Function, A::AbstractMatrix, B::AbstractMatrix,
                             alpha, beta) = unwrap_matmatmul!(C, fA, fB, A, B, alpha, beta)

@inline function unwrap_matmatmul!(C::AbstractVecOrMat, fA, fB, A::AbstractVecOrMat, B::AbstractVecOrMat,
                             alpha, beta)

    mA_axis, nA_axis = axes(fA(A))
    mB_axis, nB_axis = axes(fB(B))

    if nA_axis != mB_axis
        throw(DimensionMismatch("mul! can't contract axis $(UnitRange(nA_axis)) from A with $(UnitRange(mB_axis)) from B"))
    elseif mA_axis != axes(C,1)
        throw(DimensionMismatch("mul! got axes(C,1) == $(UnitRange(axes(C,1))), expected $(UnitRange(mA_axis)) from A"))
    elseif nB_axis != axes(C,2)
        throw(DimensionMismatch("mul! got axes(C,2) == $(UnitRange(axes(C,2))), expected $(UnitRange(nB_axis)) from B"))
    end

    # Must be sure `no_offset_view(C)` won't match signature above! 
    mul!(no_offset_view(C), fA(no_offset_view(A)), fB(no_offset_view(B)), alpha, beta)
    C
end

no_offset_view(A::Adjoint) = adjoint(no_offset_view(parent(A)))
no_offset_view(A::Transpose) = transpose(no_offset_view(parent(A)))
no_offset_view(D::Diagonal) = Diagonal(no_offset_view(parent(D)))

