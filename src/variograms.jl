using Statistics
using StatsBase
using LinearAlgebra
using Optim

# relevant links
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Geostatistical-Interpolation/Estimation-of-the-Semivariogram/index.html
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Geostatistical-Interpolation/Modeling-the-Semivariogram/index.html


function semivariogram(Z::TimeSeriesTools.AbstractRegularTimeSeries; lag_ratio=0.5, lag_max=100.0)
    Nlags = min(round(Int, lag_ratio*length(Z)*Z.Δt), round(Int, lag_max/Z.Δt))
    Npoints = length(Z)

    h = Z.Δt:Z.Δt:(Z.Δt * Nlags)
    γ = zeros(Nlags)

    # NOTE: for now, LoopVectorization.jl only supports rectangular loop indices... may be able to achieve additional
    # speedups once triangular looping is in place

    Threads.@threads for i ∈ 1:Nlags
        @inbounds for j ∈ (i+1):Npoints
            γ[i] += (Z.z[j] - Z.z[j-i])^2
        end
        @inbounds γ[i] = γ[i] / (2*(Npoints-i))
    end

    return γ, h
end


abstract type ModelVariogram end

struct SphericalVariogram{T<:Real} <: ModelVariogram
    b::T  # nugget
    C₀::T  # sill
    r::T  # range
end

function (γ::SphericalVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀ * (1.5*h/γ.r - 0.5*h^3/(γ.r^3))
    else
        return γ.b + γ.C₀
    end
end

nugget(γ::SphericalVariogram) = γ.b
sill(γ::SphericalVariogram) = γ.C₀
range(γ::SphericalVariogram) = γ.r

# this one will be much harder!
# function semivariogram(Z::AbstractGenericTimeSeries; min_window_points=10)
#     return 0
# end

# struct ModelledSemiVariogram{T <: Real}
#     nugget::T,
#     sill::T,
#     range::T,
#     method::Symbol,  #
#     θ::Vector{T}
# end



# function semivariogram_fit(h, γ)
# end
