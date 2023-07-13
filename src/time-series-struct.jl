using Unitful
using Dates
using TimeZones
#

abstract type AbstractTimeSeries end

abstract type AbstractGenericTimeSeries <: AbstractTimeSeries end
abstract type AbstractRegularTimeSeries <: AbstractTimeSeries end


struct GenericTimeSeries{T<:Real, T2, UZ<:Unitful.Units, UT<:Unitful.Units} <: AbstractGenericTimeSeries
    z::Vector{T}
    t::T2
    z_units::UZ
    t_units::UT
    start_time::ZonedDateTime
end

struct RegularTimeSeries{T<:Real, UZ<:Unitful.Units, UT<:Unitful.Units} <: AbstractRegularTimeSeries
    z::Vector{T}
    Δt::T
    z_units::UZ
    t_units::UT
    start_time::ZonedDateTime
end


"""
    times(series::AbstractGenericTimeSeries)

Given a timeseries `series` return a vector of times
"""
function times(series::AbstractGenericTimeSeries)
    return series.t
end


function times(series::AbstractRegularTimeSeries)
    return range(0.0, step=series.Δt, length=length(series.z))
end


Base.length(Z::TimeSeriesTools.AbstractTimeSeries) = length(Z.z)


# create function to convert an UncertainTimeSeries into a MultiplyUncertainTimeSeries be using a sliding window to evaluate the representativeness and variogram uncertainties
