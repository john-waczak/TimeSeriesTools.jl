module TimeSeriesTools

# Write your package code here.
include("time-series-struct.jl")
include("variograms.jl")
#include("plot_recipes.jl")
#inlcude("autocorrelation.jl")
#include("periodograms.jl")
#include("scaleograms.jl")
#inlcude("time_delay_embeddings.jl")  #<-- this could return a matrix of views of Z.z for DMD stuff
#include("koopman.jl")
#include("forecasting.jl")  # <-- arima/sarima/etc...

# time-series-struct.jl
export AbstractGenericTimeSeries, AbstractRegularTimeSeries, GenericTimeSeries, RegularTimeSeries
export UncertainGenericTimeSeries, DoublyUncertainGenericTimeSeries, TriplyUncertainGenericTimeSeries
export UncertainRegularTimeSeries, DoublyUncertainRegularTImeSeries, TriplyUncertainRegularTimeSeries
export times

# variograms.jl
export semivariogram, SphericalVariogram
export nugget, sill, γ_range
export fit_spherical_γ

end
