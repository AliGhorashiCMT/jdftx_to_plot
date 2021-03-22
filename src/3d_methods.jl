function im_epsilon_3d(lattice_vectors::Array{<:Array{<:Real, 1}, 1}, Pwannier::Array{Float64, 4}, cell_map::Array{Float64, 2} ;mesh::Int=10, spin::Int=2, histogram_width::Real = 100, )

    prefactor = e²ϵ*4*π^2*spin/3*ħ^2/(mₑ^2) #Factor of three takes into account isotropy. Spin is conventionally taken to be 2

    Epsilons = zeros(100*histogram_width)

    Vol = unit_cell_volume(lattice_vectors)

    #=  
    Note that the momentum matrix elements we'll be looking at will be in units of ev*s/angstrom.

    Therefore p^2 has units of ev^2*s^2/angstrom^2

    The integral thus has units of 1/angstrom^3 * 1/eV * ev^2 * s^2 * 1/angstrom^2 * 1/ev^2
    which is 1/angstrom^5 * ev^-1 *s^2

    The prefactor has units of ev * angstrom* 1/ev^2 * angstrom^4 * 1/s^4 * ev^2 * s^2 = ev * angstrom^5 / s^2

    Multiplying the two gives a unitless quanitity, as expected

    Where the first term is p^2, the second term is from the delta function in energy, the third term is from dividing by the Volume 
    of the unit cell, 
    =#

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            for z_mesh in 1:mesh 
                abs(momentum_matrix_elements(Pwannier, cell_map, [xmesh/mesh, y_mesh/mesh, z_mesh/mesh]))^2
            end
        end
    end
end


