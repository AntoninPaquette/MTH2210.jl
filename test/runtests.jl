using MTH2210, LinearAlgebra, Statistics, SparseArrays, Test

include("test_utils.jl")

@testset "MTH2210" begin
    
    include("test_interpolation.jl")
    include("test_non_linear.jl")
    include("test_ode.jl")

end

