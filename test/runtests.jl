using TimeSeriesTools
using Test

using Unitful
using Dates
using TimeZones

@testset "time-series-struct.jl" begin
    Δt = 0.1
    ts = collect(range(0.0, stop=60.0, step=Δt))
    zs = sin.(ts)
    z_units = u"V"  # volts
    t_units = u"s"  # seconds
    t_start = ZonedDateTime(2023, 6, 3, 12, 0, 0, tz"UTC")
    Δzs = 0.1 .* zs

    gts = GenericTimeSeries(
        zs,
        ts,
        z_units,
        t_units,
        t_start
    )

    ugts = UncertainGenericTimeSeries(
        zs,
        ts,
        z_units,
        t_units,
        t_start,
        Δzs
    )


    rts = RegularTimeSeries(
        zs,
        Δt,
        z_units,
        t_units,
        t_start
    )

    urts = UncertainRegularTimeSeries(
        zs,
        Δt,
        z_units,
        t_units,
        t_start,
        Δzs
    )

    # test that the times function works
    @test all(times(gts) .== ts)
    @test all(times(ugts) .== ts)
    @test all(times(rts) .== ts)
    @test all(times(urts) .== ts)


end



# @testset "TimeSeriesTools.jl" begin
#     # Write your tests here.
# end



