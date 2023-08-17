using TimeSeriesTools
using Statistics
using StatsBase
using LinearAlgebra
using BenchmarkTools

using CairoMakie
using CSV, DataFrames

using Unitful
using Dates, TimeZones

# using Optim,ParameterHandling
using LeastSquaresOptim,ParameterHandling



# include plotting recipe
include("makie-themes.jl")

set_theme!(mints_theme)


df_test = CSV.read(download("https://ncsa.osn.xsede.org/ees230012-bucket01/AirQualityNetwork/data/raw/Central_Hub_1/2023/03/04/MINTS_001e06318c91_IPS7100_2023_03_04.csv"), DataFrame)


typeof(df_test.pm0_1) <: AbstractVector{Float64}

Z = RegularTimeSeries(
    df_test.pm2_5,
    1.0,
    u"μg/m^3",
    u"s",
    ZonedDateTime(2023, 5, 2, 0, 0, 0, 154, tz"UTC")
)

Z2 = RegularTimeSeries(
    df_test.pm10_0,
    1.0,
    u"μg/m^3",
    u"s",
    ZonedDateTime(2023, 5, 2, 0, 0, 0, 154, tz"UTC")
)



df_test = nothing
GC.gc()



# NOTE: loopvectorization currently doesn't work. Could update in the future when they enable triangular / unstructured
# indices on nested loops.

@benchmark γ, h = semivariogram(Z, lag_ratio=0.001, lag_max=10*60)


γ, h = semivariogram(Z; lag_max=20*60)
γ2, h2 = semivariogram(Z2; lag_max=20*60)




# visualize the variogram
f = Figure()
ax = Axis(
    f[1,1],
    xlabel= "Δt (minutes)",
    ylabel = "γ(Δt)",
    title = "IPS7100 semivariogram for $(Date(Z.start_time))"
)

s1 = scatter!(
    ax,
    h ./ 60.0,
    γ,
)

s2 = scatter!(
    h2 ./ 60.0,
    γ2,
)

L = axislegend(
    ax,
    [s1, s2],
    ["PM 2.5", "PM 10.0"];
    position=:rc
)

display(f)
save("./demo_γ.pdf", f)



#------------------------
#  perform fit on γ
#------------------------

γ_params = get_reasonable_params(γ,h)

θ₀, unflatten = ParameterHandling.value_flatten(γ_params)
θ = unflatten(θ₀)
γ̂ = SphericalVariogram(θ...)

idx_fit = (h ./ 60.0) .≤ 10.0
γ_fit_spherical = fit_γ(h[idx_fit], γ[idx_fit]; method=:spherical)
γ_fit_exponential = fit_γ(h[idx_fit], γ[idx_fit]; method=:exponential)
γ_fit_gaussian = fit_γ(h[idx_fit], γ[idx_fit]; method=:gaussian)
γ_fit_circular = fit_γ(h[idx_fit], γ[idx_fit]; method=:circular)
γ_fit_cubic = fit_γ(h[idx_fit], γ[idx_fit]; method=:cubic)
γ_fit_linear = fit_γ(h[idx_fit], γ[idx_fit]; method=:linear)
γ_fit_pentaspherical = fit_γ(h[idx_fit], γ[idx_fit]; method=:pentaspherical)
γ_fit_sinehole = fit_γ(h[idx_fit], γ[idx_fit]; method=:sinehole)



fig = Figure()
ax = Axis(
    fig[1,1],
    xlabel = "Δt (minutes)",
    ylabel = "γ(Δt) for PM 2.5",
    title = "Variogram fit for PM 2.5"
)

s1 = scatter!(
    ax,
    h[idx_fit] ./ 60.0,
    γ[idx_fit],
    color=(:gray, 0.75),
    markersize = 8,
)

p1 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_spherical.(h[idx_fit]),
    linewidth=3,
)

p2 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_exponential.(h[idx_fit]),
    linewidth=3,
)

p3 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_gaussian.(h[idx_fit]),
    linewidth=3,
)

p4 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_circular.(h[idx_fit]),
    linewidth=3,
)

p5 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_cubic.(h[idx_fit]),
    linewidth=3,
)

p6 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_linear.(h[idx_fit]),
    linewidth=3,
)

p7 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_pentaspherical.(h[idx_fit]),
    linewidth=3,
)

p8 = lines!(
    h[idx_fit] ./ 60.0,
    γ_fit_sinehole.(h[idx_fit]),
    linewidth=3,
)



L = axislegend(
    ax,
    [s1, p1, p2, p3, p4, p5, p6, p7, p8],
    ["empirical", "spherical model", "exponential model", "gaussian model", "circular model", "cubic model", "linear model", "pentaspherical model", "sine hole model"];
    position=:rc
)

fig

save("./γ-fit-comparison.pdf", fig)

