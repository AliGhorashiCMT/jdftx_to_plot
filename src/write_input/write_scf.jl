function write_lattice(lattice_vectors::lattice, filename::String)
    open(filename, create=true, write=true) do io
        write(io, "lattice \\ \n")
        for lattice_row in eachrow(lattice_vectors.rvectors)
            for coord in lattice_row
                write(io, string(coord), " ")
            end
            write(io, "\\ \n" )
        end
    end
end

function write_ionpos(ionpos::ionpos, filename::String)
    open(filename, create=true, write=true) do io
        for ion in ionpos.ionpos
            for coord in ion
                write(io, string(coord))
                write(io, "  ")
            end
            write(io, "\n")
        end
    end
end

function write_scf(scf::self_consistent_field, filename::String, ionpos_filename::String, lattice_filename::String)
    open(filename, create=true, write=true, append=false) do io
        write(io, "coulomb-interaction $(scf.coulomb_interaction) \n" )
        write(io, "include  $(ionpos_filename) \n")
        write(io, "include  $(lattice_filename) \n")
        write(io, "ion-species $(scf.pseudopotential)\n")
        write(io, "elec-cutoff 20 100\n")
        write(io, "elec-initial-charge $(scf.charge)\n")
        if scf.spintype != "no-spin"
            write(io, "elec-initial-magnetization $(scf.magnetization) no \n")
        end
        write(io, "spintype $(scf.spintype)\n")
        write(io, "electronic-SCF\n")
        write(io, "dump-name $(string(filename, ".", "\$", "VAR"))\n")
        write(io, "dump End $(scf.dump)\n")
        write(io, "kpoint-folding $(scf.kpoints[1])  $(scf.kpoints[2])  $(scf.kpoints[3])\n")
        write(io, "elec-smearing Fermi $(scf.smearing)\n")

        write(io, "elec-ex-corr $(scf.xc)", "\n")
    end
end

function write_nscf(nscf::non_self_consistent_field, filename::String)
    open(filename, create=true, write=true) do io
        write(io, nscf.kpoints)
    end
end

function make_wannier_centers(scf::self_consistent_field; perturbation=10, norbitals=1)
    
    centers = []
    
    ionpos_scf = scf.ionpos

    for ion in ionpos_scf.ionpos
        coords = ion[3:5]
        for i in 1:norbitals
            push!(centers, coords + rand(Float64, 3)/perturbation)
        end
    end

    return centers

end

function write_wannier(wannier::wannier_interpolation, filename::String, scf_filename::String)

    open(filename, create=true, write=true, append=false) do io
        write(io, "include $(scf_filename) \n")
        write(io, "wannier\\ \n")
        write(io, "innerWindow $(wannier.innerWindow[1])   $(wannier.innerWindow[2])\\ \n  ")
        write(io, "outerWindow $(wannier.outerWindow[1])   $(wannier.outerWindow[2])\\ \n  ")
        write(io, "saveWfnsRealSpace", "$(wannier.saveWFNs==true ? "   yes" : "   no")", "\n")
        if wannier.phonon==true
            write(io, "\\ \n")
            for phonon_mesh in wannier.phononSupercell
                write(io, "$(phonon_mesh)  ")
            end
            write(io, "\n")
        end

        write(io, "wannier-initial-state   ", "$(string(scf_filename, ".", "\$", "VAR"))", "\n")
        write(io, "wannier-dump-name   ", "$(string(filename, ".", "\$", "VAR"))", "\n" )

        for center in wannier.wannier_centers
            write(io, "wannier-center Gaussian   ")
            for coord in center
                write(io, "$(coord)", "   ")
            end
            write(io, "\n")
        end
        write(io, "wannier-minimize niterations  $(wannier.wannier_minimize)")
    end
end


function write_phonon(phonon_params::phonon, filename::String, scf_filename::String)
    open(filename, create=true, write=true, append=false) do io
        write(io, "include  $(scf_filename)\n")
        write(io, "initial-state   ", "$(string(scf_filename, ".", "\$", "VAR"))", "\n")

        write(io, "dump-only \n\n")

        write(io, "phonon supercell  ", "$(phonon_params.supercell[1])  $(phonon_params.supercell[2])  $(phonon_params.supercell[3]) ")

    end
end