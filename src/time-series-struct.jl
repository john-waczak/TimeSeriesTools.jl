using Unitful
using Dates
using TimeZones
#

abstract type AbstractTimeSeries end

abstract type AbstractGenericTimeSeries <: AbstractTimeSeries end
abstract type AbstractRegularTimeSeries <: AbstractTimeSeries end


struct GenericTimeSeries{T<:Real, UZ<:Unitful.Units, UT<:Unitful.Units} <: AbstractGenericTimeSeries
    z::Vector{T}
    t::Vector{T}
    z_units::UZ
    t_units::UT
    start_time::ZonedDateTime
end

struct UncertainGenericTimeSeries{T<:Real, UZ<:Unitful.Units, UT<:Unitful.Units} <: AbstractGenericTimeSeries
    z::Vector{T}
    t::Vector{T}
    z_units::UZ
    t_units::UT
    start_time::ZonedDateTime
    Δz::Vector{T}
end

struct RegularTimeSeries{T<:Real, UZ<:Unitful.Units, UT<:Unitful.Units} <: AbstractRegularTimeSeries
    z::Vector{T}
    Δt::T
    z_units::UZ
    t_units::UT
    start_time::ZonedDateTime
end

struct UncertainRegularTimeSeries{T<:Real, UZ<:Unitful.Units, UT<:Unitful.Units} <: AbstractRegularTimeSeries
    z::Vector{T}
    Δt::T
    z_units::UZ
    t_units::UT
    start_time::ZonedDateTime
    Δz::Vector{T}
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



# have 2 constructors for handling units, 1 where we pass a unitful string and the second where
# where we pass a regular string and convert to unit representation.

# if no timezone provided, assume UTC

# create overloads for regular time series that will construct "t" array on-the-fly


test_unit = u"kg"
typeof(test_unit) <: Unitful.Units
