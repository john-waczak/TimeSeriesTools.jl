using TimeSeriesTools
using Statistics
using StatsBase
using LinearAlgebra
using BenchmarkTools
using Plots
using CSV, DataFrames
using Unitful
using Dates, TimeZones

using Optim, ParameterHandling

Base.length(Z::TimeSeriesTools.AbstractTimeSeries) = length(Z.z)



df_test = CSV.File(download("https://ncsa.osn.xsede.org/ees230012-bucket01/AirQualityNetwork/data/central-node-8/2023/05/02/MINTS_001e06323a37_IPS7100_2023_05_02.csv")) |> DataFrame

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

scatter(
    h ./ 60.0,
    γ,
    xlabel="Δt (minutes)",
    ylabel="γ(Δt)",
    title="IPS7100 semivariogram for $(Date(Z.start_time))",
    label="PM 2.5",
    ms=2,
    msw=0,
)

scatter!(
    h2 ./ 60.0,
    γ2,
    label="PM 10.0",
    ms=2,
    msw=0
)


savefig("./demo_γ.png")


# set up initial parameters
γ_params = (
    nugget=positive(1.5),
    sill=positive(100.0),
    range=positive(100.0)
)

γ_fit = fit_spherical_γ(h, γ, γ_params; show_trace=true)

scatter(
    h ./ 60.0,
    γ,
    xlabel="Δt (minutes)",
    ylabel="γ(Δt) for PM 2.5",
    title="IPS7100 semivariogram from $(Date(Z.start_time))",
    label="Empirical",
    ms=2,
    m=:cross,
    msw=0,
)


h_fit = 0.0:1.0:h[end]
γ_out = γ_fit.(h_fit)
plot!(
    h_fit ./ 60,
    γ_out,
    label="Spherical Model",
    lw=2,
)

scatter!([0], [nugget(γ_fit)], label="nugget", color=:grey)
vline!([γ_range(γ_fit)./60.0], ls=:dash, color=:grey, label="range")
hline!([sill(γ_fit)], ls=:dash, color=:brown, label="sill")

savefig("./demo_γ--with_fit.png")
savefig("./demo_γ--with_fit.pdf")
