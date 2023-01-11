using MTH2210, LinearAlgebra, Statistics, SparseArrays, Test

include("test_utils.jl")

# Script vérifiant les algorithmes d'interpolation
include("test_interpolation.jl")

# Script vérifiant les algorithmes de résolution d'équations non-linéaires
# include("analyse_conv_nl.jl")

# Script vérifiant les algorithmes de résolution d'équations différentielles partielles
include("test_edo.jl")
