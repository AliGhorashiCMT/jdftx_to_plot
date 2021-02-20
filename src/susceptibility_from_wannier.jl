using Plots
using PyCall
using LinearAlgebra 
using Distances

"returns the imaginary value of the polarization at frequency omega (eV) and wavevector q (inverse angstrom)"
function im_polarization(wannier_file::String, cell_map_file::String, lattice_vectors::Array{Array{Q, 1},1}, q::Array{T, 1}, μ::S; spin=2, mesh=100, histogram_width=100) where {T<:Number, Q<:Number, S<:Number}
    
    Polarization_Array=zeros(histogram_width*100)

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    np=pyimport("numpy")
    qnormalized = normalize_kvector(lattice_vectors, q)

    for i in 1:mesh
        for j in 1:mesh
            kvector=[i/mesh, j/mesh, 0]
            E1=wannier_bands(wannier_file, cell_map_file, kvector  )
            E2=wannier_bands(wannier_file, cell_map_file, kvector+qnormalized  )
            
            f1=np.heaviside( μ-E1, 0.5)
            f2=np.heaviside( μ-E2, 0.5)

            DeltaE=E2-E1
            if DeltaE>0
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*(1/mesh)^2*histogram_width
            end
        end
    end

    return Polarization_Array
end


"Applies the kramers-kronig relations onto a 1 dimensional array of numbers consisting of the imaginary value of the polarization to return the real value of polarization"
function kramers_kronig(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    #pyintegrate=pyimport("scipy.integrate")

    #pyintegrate.quad(sin, 0, max_energy)
    #ω_array=collect(1/histogram_width:1/histogram_width:max_energy)
    #sum(im_pol.*1 ./ ( ω_array.-ω.+0.001*im))*1/histogram_width
    omegaprime=collect(1:histogram_width*max_energy).*1/histogram_width
    sum(1/histogram_width*2/π*im_pol.*omegaprime./(omegaprime.^2 .- (ω+ω*0.03im)^2))

end

"returns the non-local, non-static dielectric function"
function return_2d_epsilon(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    return kramers_kronig(ω, im_pol, max_energy, histogram_width)
end