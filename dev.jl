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


df_test = CSV.File(download("https://ncsa.osn.xsede.org/ees230012-bucket01/AirQualityNetwork/data/raw/Central_Hub_1/2023/03/04/MINTS_001e06318c91_IPS7100_2023_03_04.csv")) |> DataFrame


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
γ_fit = fit_γ(h[idx_fit], γ[idx_fit])



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
    γ_fit.(h[idx_fit]),
    linewidth=3,
)

nug = scatter!(ax, [0], [nugget(γ_fit)], label="nugget", color=:grey)
rng = vlines!(ax, [γ_range(γ_fit)./60.0], ls=:dash, color=:grey, label="range")
sil = hlines!(ax, [sill(γ_fit)], ls=:dash, color=:brown, label="sill")

L = axislegend(
    ax,
    [s1, p1, nug, rng, sil],
    ["empirical", "spherical model", "nugget", "range", "sill"];
    position=:rc
)

display(fig)

save("./demo_γ--with_fit.pdf", fig)



# now let's do a demo for the other models



# # --------------------
# #    HAVOK
# # --------------------


# # Let's come up with code for generating a Hankel matrix from a time series Z.
# H = TimeDelayEmbedding(Z)

# U, σ, V = svd(H)


# # visualize the time-delay embedding
# size(V)
# scatter(V[:,1], V[:,1], V[:,3], xlabel="v₁", ylabel="v₂", zlabel="v₃",
#         ms=2,
#         msw=0,
#         msa=0,
#         marker_z = collect(1:size(V,1)) ./ (60.0)^2,
#         colormap=:viridis,
#         colorbar_title=" \n\ntime [hours]",
#         margins=5*Plots.mm,
#         label="",
#         )


# σ  # first value is *most important*

# size(U)
# size(σ)
# size(V)
# all(H .≈ U*Diagonal(σ)*V')


# plot(σ./sum(σ), xlabel="index, i", ylabel="normalized singular value σᵢ/(Σᵢσᵢ)")

# function r_cutoff(σ; ratio=0.01, rmax=15)
#     return min(sum(σ ./ sum(σ) .> ratio) + 1, rmax)
# end

# r = r_cutoff(σ)

# # truncate the matrices
# Vr = @view V[:,1:r]
# Ur = @view U[:,1:r]
# σr = @view σ[1:r]


# # compute derivatives of V matrix via fourth order central difference, ignoring last column
# dt = Z.Δt
# size(Vr)

# dVr = zeros(size(Vr,1)-5, r-1)
# for k ∈ 1:r-1, i ∈ 3:size(Vr,1)-3
#     dVr[i-2,k] = (1/(12*dt)) * (-Vr[i+2, k] + 8*Vr[i+1, k] - 8*Vr[i-1,k] + Vr[i-2,k])
# end

# @assert size(dVr,2) == r-1


# # chop off edges to size of data matches size of derivative
# X = @view Vr[3:end-3, :]
# dX = @view dVr[:,:]

# @assert size(X,2) == size(dX,2)  + 1


# # thinking of data as rows:
# # now we solve regression problem given by dX = XΞ (where Ξ is composed of our A and B matrixs)
# # solve via ordinary least squares (i.e. via normal equations)

# Ξ = (X\dX)'  # now Ξx = dx for a single column vector view

# A = Ξ[:, 1:end-1]   # State matrix A
# B = Ξ[:, end]'      # Control matrix B

# heatmap(A, yflip=true)

# eigenvecs, eigenvals = eigen(A)  # these are the modes

# plot(Vr[:,end].^2)

# size(Ξ)
# size(X)

# plot(times(Z), Z.z)
# plot!(twinx(), times(Z)[1:size(Vr,1), Vr[:,r]])
# size(Vr, 1)

# length(Z)-size(Vr,1)


# # new that we have A and B, our model should be predictive
# # starting with some initial value for v, and the control values
# # we can step forward in time:
# v₀ = Vr[1,1:r-1]
# U = Vr[:,r]




# Vr_out = 
# v_next = A*v₀ + B' .* U[1]


# p1 = plot(Vr[:,1], ylabel="v₁", label="")
# p2 = plot(Vr[:,r], ylabel="vᵣ", label="")

# p3 = plot(p1, p2, layout=(2,1))


