using Plots
using PyCall

function cell_vectors(lattice_file::String)
    run(`cat $lattice_file`);
    run(`pwd`)
end

function ion_positions(ionpos_file::String)
    run(`cat $ionpos_file`);
    run(`pwd`)
end

function plot_lattice(lattice_file::String)
    ion_position_vectors=String[]
    open(lattice_file, "r") do io
        readline(io)
        ion_position_vectors=readlines(io);
    end
end

function im_polarization()
    plot([1, 2, 3], [3, 2, 1])
end

function plot_bands(band_file::String, num_bands::Int, num_points::Int)

    reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
    exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*27.2;
    exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*27.2;

    plot(exactenergiesdown, color="black", label="", linewidth=2)
    plot!(exactenergiesup, color="purple", label="", linewidth=2)

end

function plot_phonons(cell_map::String, phononOmegaSq::String, kpoints::String)
    np=pyimport("numpy")
    cellMapPh = np.loadtxt(cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))

    eV=1/27.2
    kpointsIn = np.loadtxt(kpoints, skiprows=2, usecols=(1,2,3))
    nKin = size(kpointsIn)[1]
    #--- Fourier transform from real to k space:
    forceMatrixTilde = np.tensordot(np.exp((2im*np.pi)*np.dot(kpointsIn,transpose(cellMapPh))), forceMatrixPh, axes=1)
    #--- Diagonalize:
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTilde)
    plot(title="Phonon Dispersion", titlefontsize=20, ytickfontsize=15, sqrt.(abs.(omegaSq))/eV, linewidth=2, color="orange", legend=false, size=(800, 1000), xticks=[])

end

function pyversion()
    sys=pyimport("sys")
    print("You are currently running this version of python: $(sys.executable)")
end