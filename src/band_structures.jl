using Plots
using PyCall
using LinearAlgebra 
using Distances

#define constants for later use
eV = 1/27.2114 #in Hartrees

function plot_bands(band_file::String, num_bands::Int, num_points::Int; kwargs...)

    reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
    exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*1/eV;
    exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*1/eV;

    plot(exactenergiesdown, color="black", label="", linewidth=2; kwargs...)
    plot!(exactenergiesup, color="purple", label="", linewidth=2; kwargs...)

end

function wannier_bands(wannier_file::String, cell_map_file::String, k::Array{T, 1}) where T<:Number
    np=pyimport("numpy")
    cell_map=np.loadtxt(cell_map_file)
    Hwannier=permutedims(reshape(np.loadtxt(wannier_file), (43, 1, 1)), [1, 3, 2])
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E[1]/eV 
end