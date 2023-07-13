module TimeSeriesTools

# Write your package code here.
include("time-series-struct.jl")
include("variograms.jl")
include("embeddings.jl")
#include("plot_recipes.jl")
#inlcude("autocorrelation.jl")
#include("periodograms.jl")
#include("scaleograms.jl")
#include("koopman.jl")
#include("forecasting.jl")  # <-- arima/sarima/etc...

# time-series-struct.jl
export AbstractGenericTimeSeries, AbstractRegularTimeSeries, GenericTimeSeries, RegularTimeSeries
export UncertainGenericTimeSeries, DoublyUncertainGenericTimeSeries, TriplyUncertainGenericTimeSeries
export UncertainRegularTimeSeries, DoublyUncertainRegularTImeSeries, TriplyUncertainRegularTimeSeries
export times

# variograms.jl
export semivariogram, get_reasonable_params
export SphericalVariogram, ExponentialVariogram, GaussianVariogram, CircularVariogram, CubicVariogram, LinearVariogram, PentasphericalVariogram, SineHoleVariogram
export nugget, partialsill, sill, γ_range
export fit_spherical_γ, fit_γ

# embeddings.jl
export TimeDelayEmbedding

end
