using Statistics
using StatsBase
using LinearAlgebra
using Optim
using ParameterHandling

# relevant links
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Geostatistical-Interpolation/Estimation-of-the-Semivariogram/index.html
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Geostatistical-Interpolation/Modeling-the-Semivariogram/index.html



"""
    semivariogram(Z::TimeSeriesTools.AbstractRegularTimeSeries; lag_ratio=0.5, lag_max=100.0)

Provided a regular time series `Z`, compute the empirical [semivariogram](https://en.wikipedia.org/wiki/Variogram) γ for time lags ranging from the minimum time step Δt to a maximum value. The maximum time lag is taken to be the smaller of `lag_ratio` times the full duration of the time series and `lag_max`, a user defined maximum lag time in the time units of `Z`.

## Output

Returns a tuple  `(γ, h)` containing the (scaled) variances and time lags.
"""
function semivariogram(Z::TimeSeriesTools.AbstractRegularTimeSeries; lag_ratio=0.5, lag_max=100.0)
    # sanity check for supplied arguments
    if lag_ratio ≤ 0 || 1 ≤ lag_ratio
        throw(DomainError(lag_ratio, "argument ∉ (0,1)"))
    end

    if lag_max ≤ 0
        throw(DomainError(lag_max, "argument must be positive"))
    end


    Nlags = min(round(Int, lag_ratio*length(Z)*Z.Δt), round(Int, lag_max/Z.Δt))

    Npoints = length(Z)

    h = Z.Δt:Z.Δt:(Z.Δt * Nlags)
    γ = zeros(Nlags)

    # NOTE: for now, LoopVectorization.jl only supports rectangular loop indices... may be able to achieve additional
    # speedups once triangular looping is in place


    # add a flag to choose whether or not to use multithreading
    Threads.@threads for i ∈ 1:Nlags
        @inbounds for j ∈ (i+1):Npoints
            γ[i] += (Z.z[j] - Z.z[j-i])^2
        end
        @inbounds γ[i] = γ[i] / (2*(Npoints-i))
    end

    return γ, h
end




abstract type ModelVariogram end

struct SphericalVariogram{T} <: ModelVariogram
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
γ_range(γ::SphericalVariogram) = γ.r


"""
    fit_spherical_γ(h, γ, params; optargs...)

Given an empirical semi-variogram γ with lags h, fit a spherical variogram model to the data with parameters `params`. The arguments `optargs` are passed to the `optimize` function.
"""
function fit_spherical_γ(h, γ, params; optargs...)
    θ₀, unflatten = ParameterHandling.value_flatten(params)

    function loss(θ)
        ps = unflatten(θ)
        γ̂ = SphericalVariogram(ps.nugget, ps.sill, ps.range).(h)
        return sqrt(mean((γ .- γ̂).^2))
    end

    opt = optimize(loss, θ₀; optargs...)

    ps_fitted = unflatten(opt.minimizer)

    return SphericalVariogram(ps_fitted.nugget, ps_fitted.sill, ps_fitted.range)
end








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
