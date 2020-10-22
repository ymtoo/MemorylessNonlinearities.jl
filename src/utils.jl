minmaxrescale(x::AbstractArray{T}, minval, maxval) where {T<:Real} = minval .+ (
    maxval - minval) .* (x .-minimum(x)) ./ (maximum(x) - minimum(x))