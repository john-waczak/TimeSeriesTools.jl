using TimeSeriesTools
using Test

using Unitful
using Dates
using TimeZones


Δt = 0.1
ts = collect(range(0.0, stop=60.0, step=Δt))
zs = sin.(ts)
z_units = u"V"  # volts
t_units = u"s"  # seconds
t_start = ZonedDateTime(2023, 6, 3, 12, 0, 0, tz"UTC")
Δzs = 0.1 .* zs

Z = RegularTimeSeries(
    zs,
    Δt,
    z_units,
    t_units,
    t_start
)



@testset "time-series-struct.jl" begin
    gts = GenericTimeSeries(
        zs,
        ts,
        z_units,
        t_units,
        t_start
    )

    rts = RegularTimeSeries(
        zs,
        Δt,
        z_units,
        t_units,
        t_start
    )

    # test that the times function works
    @test all(times(gts) .== ts)
    @test all(times(rts) .== ts)

    @test length(gts) == length(rts)

end


@testset "variograms.jl" begin

    # check that supplied arguments are reasonable
    @test_throws DomainError semivariogram(Z, lag_ratio=-1.0)
    @test_throws DomainError semivariogram(Z, lag_max=-1.0)
end



@testset "embeddings.jl" begin
    out = [
        1 2 3
        2 3 4
        3 4 5
    ]

    out2 = [
        3 4 5
        2 3 4
        1 2 3
    ]

    @test all(TimeDelayEmbedding(1:5, nrow=3) .== out)
    @test all(TimeDelayEmbedding(1:5, nrow=3, method=:backward) .== out2)

end



# @testset "TimeSeriesTools.jl" begin
#     # Write your tests here.
# end



