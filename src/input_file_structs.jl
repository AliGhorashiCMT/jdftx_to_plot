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
    wannier_centers::Array{Array{Int64, 1}, 1}
    saveWFNs::Bool
    phonon::Bool
    phononSupercell::Array{Int64, 1}
    innerWindow::Array{T, 1} where T<:Number
    outerWindow::Array{S, 1} where S<:Number
    wannier-minimize::Int
end

wannier_interpolation(wannier_centers, inner_Window, outer_Window)=wannier_interpolation(wannier_centers, false, false, [0, 0, 0], innerWindow, outerWindow, 10000)
