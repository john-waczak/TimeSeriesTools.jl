module TimeSeriesTools

# Write your package code here.
include("time-series-struct.jl")
include("variograms.jl")
# time-delay-embeddings.jl <--- good for DMD style stuff perhaps?


export GenericTimeSeries, UncertainGenericTimeSeries, RegularTimeSeries, UncertainRegularTimeSeries
export times

end
