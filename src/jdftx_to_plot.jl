module jdftx_to_plot

#= We define all constants required for future calculations.
hbar is given in ev*seconds. 
c, the speed of light, is given in angstroms/second
α is the fine structure constant (unitless)
e^2ϵ is actually e^2/ϵ and is given in units of ev*angstrom
the constant eV is the conversion from hartrees to eV.
Note that jdftx output is always in hartrees (atomic units)
=#
const ħ = 6.6e-16
const c = 3e18
const α = 1/137
const e²ϵ = 4π*ħ*c*α  
const bohrtoangstrom = 0.529177
const eV = 1/27.2114 

export ħ,c, e²ϵ, bohrtoangstrom

include("input_file_structs.jl")
export self_consistent_field, non_self_consistent_field, wannier_interpolation,
lattice, ionpos, phonon

#= 
Properties of the unit cell, the real space lattice and the reciprocal lattice. 
Methods provided for normalizing the kvectors in the basis of reciprocal lattice vectors (necessary for jdftx).
=#
include("cell_properties.jl")
export unnormalize_kvector, normalize_kvector, brillouin_zone_area,
in_wigner_seitz, in_brillouin, reciprocal_vectors, ion_positions, plot_lattice, 
cell_vectors

include("phonon_properties.jl")
export plot_phonons, phonon_dispersion

#=
A multitude of methods for calculating the density of states- either directly from DFT output or from wannier 
tight binding data. Functions also provided for the density of states of phonons. Lastly, a handy method has been provided
for calculation of the chemical potential at arbitrary filling. 
=#
include("density_of_states.jl")
export density_of_states, density_of_states_wannier, find_chemical_potential, phonon_density_of_states

#=
Methods to find the susceptibility of materials and their logitudinal dielectric response (real and imaginary).
A variety of different methods are provided for ease of cross checking results. For instance, one may opt to find the real part of the susceptibility 
either through kramers kronig of the imaginary susceptibility or directly from 2d integration.
For Kramers-Kronig, we've also provided methods to find the imaginary susceptibility from the real susceptibility. 
This may be used to ensure results are reliable. 
=#
include("susceptibility_from_wannier.jl")
export im_polarization, kramers_kronig, kramers_kronig_scipy, kramers_kronig_quadgk, im_polarization_cubature, 
return_2d_epsilon, return_2d_epsilon_scipy, direct_epsilon,
direct_epsilon_cubature, return_2d_epsilon_quadgk,
kramers_kronig_reverse_scipy, kramers_kronig_reverse_quadgk, kramers_kronig_reverse

#=
Methods to plot band structures- either from direct DFT data or from wannier tight binding data 
=#
include("band_structures.jl")
export wannier_bands, wannier_vectors, plot_bands, hwannier 

include("analytic_models.jl")

include("export_wannier_hamiltonians.jl")

#=
smoothing functions- useful for kramers kronig calculations for which a smooth imaginary susceptibility is preferable 
for reliable numerics. 
=#
include("smooth.jl")
export smooth

include("matrix_elements.jl")
export phmatrixelements

#= 
Methods to create supercells/large defect lattices using an underlying smaller unit cell
=#
include("supercell.jl")
export make_supercell, make_defectcell

#=
Methods to write DFT input files. Note that these input files are specifically written with JDFTX in mind. 
=#
include("./write_input/write_scf.jl")
export write_scf, write_nscf, write_wannier, write_ionpos, write_lattice, write_phonon

#=
Methods to calculate damping of plasmons up to second order in phonon interactions 
=#
include("./loss_calculations/plasmon_losses.jl")
export landau_damping, first_order_damping, second_order_damping


end # module
