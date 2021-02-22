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
            end
            println()
        end
    end
end

function write_scf(scf::self_consistent_field, filename::String)
    open(filename, create=true, write=true) do io
        write(io, scf.xc)
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