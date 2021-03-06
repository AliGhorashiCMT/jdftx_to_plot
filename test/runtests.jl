#using jdftx_to_plot
using Test, PyCall, jdftx_to_plot

@testset "jdftx_to_plot" begin
    include("wannier_bands_tests.jl")
    include("cellsizes.jl")
    include("dos.jl")
    include("plasmons.jl")
    include("matrix_elements_tests.jl")
end
