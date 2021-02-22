struct lattice
    ionpos::Array{Float64, 2}
end

struct ionpos
    ionpos::Array{Array{Any, 1}, 1}
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
    dump::String #Array{String, 1}
end

self_consistent_field(xc, kpoints, lattice, ionpos)=self_consistent_field(xc, kpoints, lattice, ionpos,"GBRV/\$ID_pbsesol.uspp", 0, "no-spin", 0, 0.00001, "ElecDensity" )

struct non_self_consistent_field
    scf::self_consistent_field
    kpoints::Array{Int, 1}
end

struct wannier_interpolation
    wannier_centers::Array{Int64, 2}
    phononSupercell::Array{Int64, 1}
end

