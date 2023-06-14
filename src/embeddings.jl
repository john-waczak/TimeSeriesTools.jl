

"""
    TimeDelayEmbedding(Z::AbstractVector; nrow=100, method=:Forward)

NOTE: The resulting matrix treats each column in the returned matrix `H` as an embedded featured vector.
"""
function TimeDelayEmbedding(Z::AbstractVector; nrow=100, method=:Forward)
    if nrow ≥ length(Z)
        throw(DomainError(nrow, "nrow must be less than the full length of Z"))
    end

    if method ∉ [:Forward, :Backward]
        throw(ArgumentError("method must be one of :Forward or :Backwards"))
    end

    ncol = length(Z) - nrow + 1

    # pre-allocate output
    H = zeros(nrow, ncol)

    for k ∈ 1:ncol
        H[:, k] = Z[k:k+nrow-1]
    end

    if method == :Backward
        H .= H[end:-1:1, :]
    end

    return H
end


function TimeDelayEmbedding(Z::AbstractRegularTimeSeries; kwargs...)
    return TimeDelayEmbedding(Z.z; kwargs...)
end
