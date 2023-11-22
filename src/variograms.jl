using Statistics
using StatsBase
using LinearAlgebra
#using LeastSquaresOptim
#using ParameterHandling
using Optimization
using OptimizationOptimJL



# https://github.com/gokulbalagopal/Temporal-Variograms/blob/master/firmware/non_linear_fit.jl#L45


# relevant links
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Geostatistical-Interpolation/Estimation-of-the-Semivariogram/index.html
# https://www.geo.fu-berlin.de/en/v/soga/Geodata-analysis/geostatistics/Geostatistical-Interpolation/Modeling-the-Semivariogram/index.html



# see this link for list of theoretical models:
# https://ml-gis-service.com/index.php/2022/04/01/geostatistics-theoretical-variogram-models/

# also:
# https://scikit-gstat.readthedocs.io/en/latest/reference/models.html#cubic-modeG

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



function semivariogram(Zs, Δt; lag_ratio=0.5, lag_max=100.0)
    # sanity check for supplied arguments
    if lag_ratio ≤ 0 || 1 < lag_ratio
        throw(DomainError(lag_ratio, "argument ∉ (0,1)"))
    end

    if lag_max ≤ 0
        throw(DomainError(lag_max, "argument must be positive"))
    end


    Nlags = min(round(Int, lag_ratio*length(Zs)*Δt), round(Int, lag_max/Δt))

    Npoints = length(Zs)

    h = Δt:Δt:(Δt * Nlags)
    γ = zeros(Nlags)

    # NOTE: for now, LoopVectorization.jl only supports rectangular loop indices... may be able to achieve additional
    # speedups once triangular looping is in place


    # add a flag to choose whether or not to use multithreading
    Threads.@threads for i ∈ 1:Nlags
        @inbounds for j ∈ (i+1):Npoints
            γ[i] += (Zs[j] - Zs[j-i])^2
        end
        @inbounds γ[i] = γ[i] / (2*(Npoints-i))
    end

    return γ, h
end





"""
    get_reasonable_params(γ,h)

Given an empirical semivariogram `γ` with time lags `h`, return reasonable default parameters for fitting a variogram model using the ParameterHandling format with positive constraints. Nugget is defaulted to 0.01 of maximum γ, partialsill is set to 85% of max γ value, and range is set to 10% of maximum lag.
"""
function get_reasonable_params(γ,h)
    # return (
    #     nugget = positive(0.05*maximum(γ)),
    #     partialsill = positive(0.9*maximum(γ)),
    #     range = positive(0.15*maximum(h))
    # )
    γ_min, γ_max = extrema(γ)
    Δγ = γ_max - γ_min
    return [γ_min + 0.1*Δγ, γ_min + 0.75*Δγ, 0.10*maximum(h)]
end




abstract type ModelVariogram end

struct SphericalVariogram{T} <: ModelVariogram
    b::T  # nugget
    C₀::T  # partialsill
    r::T  # range
end

function (γ::SphericalVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀ * (1.5*h/γ.r - 0.5*h^3/(γ.r^3))
    else
        return γ.b + γ.C₀
    end
end


struct ExponentialVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::ExponentialVariogram)(h)
    return γ.b + γ.C₀*(1-exp(-h/γ.r))
end


struct GaussianVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::GaussianVariogram)(h)
    return γ.b + γ.C₀*(1-exp(-h^2/(γ.r)^2))
end



struct CircularVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::CircularVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀*(1-((2/π)*acos(h/γ.r)) + (2*h/(π*γ.r)) * sqrt(1-(h/γ.r)^2))
    else
        return γ.b + γ.C₀
    end
end


struct CubicVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::CubicVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀*(7*(h/γ.r)^2 - 8.75*(h/γ.r)^3 + 3.5*(h/γ.r)^5 - 0.75*(h/γ.r)^7 )
    else
        return γ.b + γ.C₀
    end
end


struct LinearVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::LinearVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀*h/γ.r
    else
        return γ.b + γ.C₀
    end
end


struct PentasphericalVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::PentasphericalVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀*((15/8)*(h/γ.r) - (5/4)*(h/γ.r)^3 + (3/8)*(h/γ.r)^5)
    else
        return γ.b + γ.C₀
    end
end


struct SineHoleVariogram{T} <: ModelVariogram
    b::T   # nugget
    C₀::T  # sill
    r::T   # range
end

function (γ::SineHoleVariogram)(h)
    if h ≤ γ.r
        return γ.b + γ.C₀*(1 - (sin(π*h/γ.r))/(π*h/γ.r))
    else
        return γ.b + γ.C₀
    end
end



nugget(γ::ModelVariogram) = γ.b
sill(γ::ModelVariogram) = γ.C₀ + γ.b
partialsill(γ::ModelVariogram) = γ.C₀
γ_range(γ::ModelVariogram) = γ.r



method_lookup = Dict(
    :spherical => SphericalVariogram,
    :exponential => ExponentialVariogram,
    :gaussian => GaussianVariogram,
    :circular => CircularVariogram,
    :cubic => CubicVariogram,
    :linear => LinearVariogram,
    :pentaspherical => PentasphericalVariogram,
    :sinehole => SineHoleVariogram
)



function fit_γ(h, γ, u0; method=:spherical)

    f = method_lookup[method]

    function loss(u, p)
        γ̂ = f(u[1], u[2], u[3])
        return sum(x->x^2, γ̂.(h) .- γ)
    end

    optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optf, u0, nothing, lb = [0.0, minimum(γ), minimum(h)], ub = [maximum(γ), maximum(γ), maximum(h)])

    sol = solve(prob, BFGS())

    return f(sol[1], sol[2], sol[3])
end



function fit_γ(h, γ; method=:spherical)

    u0 = get_reasonable_params(γ,h)
    f = method_lookup[method]

    function loss(u, p)
        γ̂ = f(u[1], u[2], u[3])
        return sum(x->x^2, γ̂.(h) .- γ)
    end

    optf = OptimizationFunction(loss, Optimization.AutoForwardDiff())

    println("u0: $(u0)")
    println("lb: $([0.0, minimum(γ), minimum(h)])")
    println("ub: $([maximum(γ), maximum(γ), maximum(h)])")

    prob = OptimizationProblem(optf, u0, nothing, lb = [0.0, minimum(γ), minimum(h)], ub = [maximum(γ), maximum(γ), maximum(h)])

    sol = solve(prob, BFGS())

    return f(sol[1], sol[2], sol[3])
end

