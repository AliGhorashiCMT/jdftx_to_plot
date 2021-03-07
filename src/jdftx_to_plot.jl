module jdftx_to_plot
using PyCall

const ħ = 6.6e-16
const c = 3e18
const α = 1/137
const e²ϵ = 4π*ħ*c*α  
const bohrtoangstrom = 0.529177

export(ħ)
export(c)
export(e²ϵ)
export(bohrtoangstrom)

include("input_file_structs.jl")

include("cell_properties.jl")
include("phonon_properties.jl")
include("density_of_states.jl")
include("susceptibility_from_wannier.jl")
include("band_structures.jl")
include("analytic_models.jl")
include("export_wannier_hamiltonians.jl")
include("smooth.jl")
include("matrix_elements.jl")

include("supercell.jl")

export(self_consistent_field)
export(non_self_consistent_field)
export(wannier_interpolation)
export(lattice)
export(ionpos)
export(phonon)

include("./write_input/write_scf.jl")
export(write_scf)
export(write_nscf)
export(write_wannier)
export(write_ionpos)
export(write_lattice)
export(write_phonon)


include("./loss_calculations/plasmon_losses.jl")
export(landau_damping)
export(first_order_damping)
export(second_order_damping)

export(cell_vectors)
export(ion_positions)
export(plot_lattice)
export(plot_bands)
export(plot_phonons)
export(kramers_kronig)
export(wannier_bands)
export(reciprocal_vectors)
export(brillouin_zone_area)
export(im_polarization)
export(normalize_kvector)
export(density_of_states)
export(density_of_states_wannier)
export(in_wigner_seitz)
export(return_2d_epsilon)
export(kramers_kronig_scipy)
export(return_2d_epsilon_scipy)
export(direct_epsilon)
export(direct_epsilon_cubature)
export(smooth)
export(kramers_kronig_quadgk)
export(return_2d_epsilon_quadgk)
export(im_polarization_cubature)

export(supercell)

export(phmatrixelements)


export(phonon_density_of_states)
export(phonon_dispersion)

export(unnormalize_kvector)
export(wannier_vectors)
export(find_chemical_potential)

export(in_wigner_seitz)
export(in_brillouin)
export(hwannier)

greet() = print("Hello World!")




export(greet)

end # module
