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
    nugget=positive(0.01),
    sill=positive(0.1),
    range=positive(100.0)
)

γ_fit = fit_spherical_γ(h, γ, γ_params; show_trace=true)

println(γ_fit)

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
hline!([sill(γ_fit) + nugget(γ_fit)], ls=:dash, color=:brown, label="sill")


println(sqrt(nugget(γ_fit)))

mean(Z.z)

savefig("./demo_γ--with_fit.png")
savefig("./demo_γ--with_fit.pdf")



# --------------------
#    HAVOK
# --------------------


# Let's come up with code for generating a Hankel matrix from a time series Z.
typeof(Z)

function TimeDelayEmbedding(Z::AbstractRegularTimeSeries; nrow=100)
    # pre-allocate resulting matrix H
    H = zeros(nrow, length(Z)-nrow)

    # fill up H with our data
    ncol = length(Z) - nrow
    for k∈1:ncol
        H[:, k] = Z.z[k:end-ncol+k-1]
    end

    # here we think of rows of H as our time series (shifted).
    return H
end


H = TimeDelayEmbedding(Z)

U, σ, V = svd(H)


# visualize the time-delay embedding
size(V)
scatter(V[:,1], V[:,1], V[:,3], xlabel="v₁", ylabel="v₂", zlabel="v₃",
        ms=2,
        msw=0,
        msa=0,
        marker_z = collect(1:size(V,1)) ./ (60.0)^2,
        colormap=:viridis,
        colorbar_title=" \n\ntime [hours]",
        margins=5*Plots.mm,
        label="",
        )


σ  # first value is *most important*

size(U)
size(σ)
size(V)
all(H .≈ U*Diagonal(σ)*V')


plot(σ./sum(σ), xlabel="index, i", ylabel="normalized singular value σᵢ/(Σᵢσᵢ)")

function r_cutoff(σ; ratio=0.01, rmax=15)
    return min(sum(σ ./ sum(σ) .> ratio) + 1, rmax)
end

r = r_cutoff(σ)

# truncate the matrices
Vr = @view V[:,1:r]
Ur = @view U[:,1:r]
σr = @view σ[1:r]


# compute derivatives of V matrix via fourth order central difference, ignoring last column
dt = Z.Δt
size(Vr)

dVr = zeros(size(Vr,1)-5, r-1)
for k ∈ 1:r-1, i ∈ 3:size(Vr,1)-3
    dVr[i-2,k] = (1/(12*dt)) * (-Vr[i+2, k] + 8*Vr[i+1, k] - 8*Vr[i-1,k] + Vr[i-2,k])
end

@assert size(dVr,2) == r-1


# chop off edges to size of data matches size of derivative
X = @view Vr[3:end-3, :]
dX = @view dVr[:,:]

@assert size(X,2) == size(dX,2)  + 1


# thinking of data as rows:
# now we solve regression problem given by dX = XΞ (where Ξ is composed of our A and B matrixs)
# solve via ordinary least squares (i.e. via normal equations)

Ξ = (X\dX)'  # now Ξx = dx for a single column vector view

A = Ξ[:, 1:end-1]   # State matrix A
B = Ξ[:, end]'      # Control matrix B

heatmap(A, yflip=true)

eigenvecs, eigenvals = eigen(A)  # these are the modes

plot(Vr[:,end].^2)

size(Ξ)
size(X)

plot(times(Z), Z.z)
plot!(twinx(), times(Z)[1:size(Vr,1), Vr[:,r]])
size(Vr, 1)

length(Z)-size(Vr,1)


# new that we have A and B, our model should be predictive
# starting with some initial value for v, and the control values
# we can step forward in time:
v₀ = Vr[1,1:r-1]
U = Vr[:,r]




Vr_out = 
v_next = A*v₀ + B' .* U[1]


p1 = plot(Vr[:,1], ylabel="v₁", label="")
p2 = plot(Vr[:,r], ylabel="vᵣ", label="")

p3 = plot(p1, p2, layout=(2,1))


