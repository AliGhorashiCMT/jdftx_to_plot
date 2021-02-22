function write_lattice(lattice::lattice, filename::String)
    open(filename, create=true, write=true) do io
        write(io, ionpos.lattice)
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

        write(io, "include  $(ionpos_filename) \n")
        write(io, "include  $(lattice_filename) \n")
        write(io, "ion-species $(scf.pseudopotential)\n")
        write(io, "elec-cutoff 20 100\n")
        write(io, "elec-initial-charge $(scf.charge)\n")
        write(io, "spintype $(scf.spintype)\n")
        write(io, "electronic-SCF\n")
        write(io, "dump-name $(string(filename, ".", "\$", "VAR"))\n")
        write(io, "dump End $dump\n")
        write(io, "kpoint-folding $(scf.kpoints[1]), $(scf.kpoints[2]), $(scf.kpoints[3])\n")
        write(io, "elec-smearing Fermi $(scf.smearing)\n")

        write(io, scf.xc, "\n")
    end
end

function write_nscf(nscf::non_self_consistent_field, filename::String)
    open(filename, create=true, write=true) do io
        write(io, nscf.kpoints)
    end
end

function write_wannier(wannier::wannier_interpolation, filename::String)
    open(filename, create=true, write=true) do io
        write(io, wannier.wannier_centers)
    end
end