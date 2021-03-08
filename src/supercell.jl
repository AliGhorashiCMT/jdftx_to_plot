function make_supercell(small_lattice::lattice, small_ionpos::ionpos, cell_mult::Array{Int, 1})
    
    mult1, mult2, mult3 = cell_mult
    
    supercell_ionpos = []
    for ion_positions in small_ionpos.ionpos
        new_base_position = ion_positions[3:end].*(1/mult1, 1/mult2, 1/mult3, 1)
        label_1 = ion_positions[1]
        label_2 = ion_positions[2]
        for i in 0:mult1-1
            for j in 0:mult2-1
                for k in 0:mult3-1
                    new_pos=[] 
                    push!(new_pos, label_1, label_2)

                    append!(new_pos, new_base_position +[i/mult1, j/mult2, k/mult3, 0] )
                    push!(supercell_ionpos, new_pos)

                end
            end
        end
    end

    return lattice([small_lattice.rvectors[:, 1]*cell_mult[1] small_lattice.rvectors[:, 2]*cell_mult[2] small_lattice.rvectors[:, 3]*cell_mult[3]]), supercell_ionpos 
end

function make_defectcell(small_lattice::lattice, small_ionpos::ionpos, cell_mult::Array{Int, 1}, defect_atom::String)

    defect_lattice, supercell_ionpos =  make_supercell(small_lattice, small_ionpos, cell_mult)
    supercell_ionpos[1][2] = defect_atom

    return defect_lattice, supercell_ionpos

end