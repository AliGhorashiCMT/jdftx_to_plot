using Plots
using PyCall
using LinearAlgebra 
using Distances
using HCubature
using QuadGK


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

"Applies the kramers-kronig relations but with scipy's cauchy weight; kwargs for scipy.integrate.quad supported"
function kramers_kronig_scipy(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q, max_energy_integration; kwargs...) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    pyintegrate=pyimport("scipy.integrate")
    interpol=pyimport("scipy.interpolate")
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(max_energy-1/histogram_width), im_pol)
    
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)

    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]

end

function kramers_kronig_quadgk(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q, max_energy_integration; δ=.1, kwargs...) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    
    interpol=pyimport("scipy.interpolate")
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(max_energy-1/histogram_width), im_pol)
    
   inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime^2-(ω+1im*δ)^2)

    return real(quadgk(inner_function, 0, max_energy_integration; kwargs...)[1])

end


function epsilon_integrand(wannier_file, cell_map_file, k₁, k₂, q, μ, ω, ϵ; spin=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(wannier_file, cell_map_file, kvector,  )
    ϵ₂ =wannier_bands(wannier_file, cell_map_file, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    real(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

function epsilon_integrand_imaginary(wannier_file, cell_map_file, k₁, k₂, q, μ, ω, ϵ; spin=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(wannier_file, cell_map_file, kvector,  )
    ϵ₂ =wannier_bands(wannier_file, cell_map_file, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    imag(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

function direct_epsilon(wannier_file::String, cell_map_file::String, lattice_vectors::Array{Array{Q, 1},1}, q::Array{T, 1}, ω::R, μ::S; spin=1, ϵ=0.01, kwargs...) where {T<:Number, Q<:Number, S<:Number, R<:Number}
    
    kwargsdict=Dict()

    for kwarg in kwargs
        push!(kwargsdict, kwarg.first => kwarg.second)
    end

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))

    pyintegration=pyimport("scipy.integrate")

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*pyintegration.nquad((k₁, k₂) -> epsilon_integrand(wannier_file, cell_map_file, k₁, k₂, qnormalized, μ, ω, ϵ, spin=spin), [[0, 1], [0, 1]], opts=kwargsdict)[1]

    1-90.5/qabs*polarization

end

"Direct 2D integration for Epsilon with HCubature"
function direct_epsilon_cubature(wannier_file::String, cell_map_file::String, lattice_vectors::Array{Array{Q, 1},1}, q::Array{T, 1}, ω::R, μ::S; spin=1, ϵ=0.01, kwargs...) where {T<:Number, Q<:Number, S<:Number, R<:Number}

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))

    pyintegration=pyimport("scipy.integrate")

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand(wannier_file, cell_map_file, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]

    1-90.5/qabs*polarization

end

"Find the imaginary value of polarization through hcubature "
function im_polarization_cubature(wannier_file::String, cell_map_file::String, lattice_vectors::Array{Array{Q, 1},1}, q::Array{T, 1}, ω::R, μ::S; spin=1, ϵ=0.01, kwargs...) where {T<:Number, Q<:Number, S<:Number, R<:Number}

    qnormalized = normalize_kvector(lattice_vectors, q)

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand_imaginary(wannier_file, cell_map_file, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]

    return polarization

end


"returns the non-local, non-static dielectric function"
function return_2d_epsilon(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    return kramers_kronig(ω, im_pol, max_energy, histogram_width)
end

"returns the non-local, non-static dielectric function using scipy functionality"
function return_2d_epsilon_scipy(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    return kramers_kronig_scipy(ω, im_pol, max_energy, histogram_width)
end

function return_2d_epsilon_quadgk(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q; δ=0.1, kwargs... ) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    return kramers_kronig_quadgk(ω, im_pol, max_energy, histogram_width; δ, kwargs...)  
end