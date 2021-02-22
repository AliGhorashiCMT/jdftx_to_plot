struct lattice

end

struct ionpos

end

struct self_consistent_field
    xc::String
    kpoints::Array{Int, 1}
    lattice::lattice
    ionpos::ionpos
    pseudopotential::String
    charge::T where T<:Number
    spintype::String
    magnetization::R where R<:Number
    smearing::S where S<:Number
    dump::Array{String, 1}
    Fermi::Float64
end

self_consistent_field(xc, kpoints, lattice)=self_consistent_field(xc, kpoints, lattice, ionpos,"GBRV/\$ID_pbsesol.uspp" 0, "no-spin", 0)

struct non_self_consistent_field
    scf::self_consistent_field
    kpoints::Array{Int, 1}
end

struct wannier_interpolation
    wannier_centers::Array{Int64, 2}
    phononSupercell::Array{Int64, 1}
end

