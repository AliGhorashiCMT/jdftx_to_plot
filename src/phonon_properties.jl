using Plots
using PyCall
using LinearAlgebra 
using Distances

"Plots the phonon band dispersion at the kpoints supplied"
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

"Give phonon dispersion at individual kpoints"
function phonon_dispersion(cell_map::String, phononOmegaSq::String, qnorm::Array{T, 1}), 
    np=pyimport("numpy")
    cellMapPh = np.loadtxt(cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))

    eV=1/27.2
    #--- Fourier transform from real to k space:
    forceMatrixTildeq = np.tensordot(np.exp((2im*np.pi)*np.dot(qnorm,transpose(cellMapPh))), forceMatrixPh, axes=1)
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTildeq)
    return sqrt.(abs.(omegaSq))/eV
end