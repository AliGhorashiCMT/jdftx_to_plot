"""
Returns the imaginary value of the polarization at frequency omega (eV) and wavevector q (inverse angstrom).
Several methods are provided. Wannier and cell map data may be given either through file names or through passing in 
HWannier and cell-map as dim 3 and dim 2 arrays of floats, respectively.
"""
function im_polarization(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; spin::Int = 1, mesh::Int = 100, histogram_width::Real = 100) 
    
    Polarization_Array=zeros(histogram_width*100)

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
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
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
            end
        end
    end

    return Polarization_Array
end

function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; spin::Int=1, mesh::Int=100, histogram_width::Real=100) 
    
    Polarization_Array=zeros(histogram_width*100)

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)

    for i in 1:mesh
        for j in 1:mesh
            kvector=[i/mesh, j/mesh, 0]
            E1 = wannier_bands(HWannier, cell_map, kvector  )
            E2 = wannier_bands(HWannier, cell_map, kvector+qnormalized  )
            
            f1=np.heaviside( μ-E1, 0.5)
            f2=np.heaviside( μ-E2, 0.5)

            DeltaE=E2-E1
            if DeltaE>0
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
            end
        end
    end

    return Polarization_Array
end


function im_polarization(wannier_file::String, cell_map_file::String, lattvectors::lattice, q::Array{<:Real, 1}, μ::Real; spin::Int = 1, mesh::Int = 100, histogram_width::Real = 100)
    
    Polarization_Array=zeros(histogram_width*100)

    lattice_vectors = [lattvectors.rvectors[:, 1]*bohrtoangstrom, lattvectors.rvectors[:, 2]*bohrtoangstrom, lattvectors.rvectors[:, 3]*bohrtoangstrom]

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
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
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
            end
        end
    end

    return Polarization_Array
end


function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattvectors::lattice, q::Array{<:Real, 1}, μ::Real; spin::Int=1, mesh::Int=100, histogram_width::Int=100) 
    
    Polarization_Array=zeros(histogram_width*100)

    lattice_vectors = [lattvectors.rvectors[:, 1]*bohrtoangstrom, lattvectors.rvectors[:, 2]*bohrtoangstrom, lattvectors.rvectors[:, 3]*bohrtoangstrom]

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)

    for i in 1:mesh
        for j in 1:mesh
            kvector=[i/mesh, j/mesh, 0]
            E1=wannier_bands(HWannier, cell_map, kvector  )
            E2=wannier_bands(Hwannier, cell_map, kvector+qnormalized  )
            
            f1=np.heaviside( μ-E1, 0.5)
            f2=np.heaviside( μ-E2, 0.5)

            DeltaE=E2-E1
            if DeltaE>0
                Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*(1/mesh)^2*histogram_width*spin
            end
        end
    end

    return Polarization_Array
end


function im_polarization(wannier_file::String, cell_map_file::String, nbands::Int, valence_bands::Int, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; spin::Int=1, mesh::Int=100, histogram_width::Int=100) 
    
    Polarization_Array=zeros(histogram_width*100)

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)

    for i in 1:mesh
        for j in 1:mesh
            kvector=[i/mesh, j/mesh, 0]
            E1=wannier_bands(wannier_file, cell_map_file, kvector, nbands  )
            E2=wannier_bands(wannier_file, cell_map_file, kvector+qnormalized, nbands  )
            
            V1=wannier_vectors(wannier_file, cell_map_file, kvector, nbands  )
            V2=wannier_vectors(wannier_file, cell_map_file, kvector+qnormalized, nbands  )

            for lower in 1:valence_bands+1
                for upper in valence_bands+1:nbands
                    Elower = E1[lower]
                    Eupper = E2[upper]
                    overlap=(np.abs(np.dot(V1[:, lower], np.conj(V2[:, upper]))))^2;
                    f1=np.heaviside( μ-Elower, 0.5)
                    f2=np.heaviside( μ-Eupper, 0.5)

                    DeltaE=Eupper-Elower
                    if DeltaE>0
                        Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*overlap*(1/mesh)^2*histogram_width*spin
                    end
                end
            end
        end
    end

    return Polarization_Array
end


function im_polarization(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Int, valence_bands::Int, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; spin::Int=1, mesh::Int=100, histogram_width::Int=100) 
    
    Polarization_Array=zeros(histogram_width*100)

    V=(2π)^2/brillouin_zone_area(lattice_vectors)
    qnormalized = normalize_kvector(lattice_vectors, q)

    for i in 1:mesh
        for j in 1:mesh
            kvector=[i/mesh, j/mesh, 0]
            E1=wannier_bands(HWannier, cell_map, kvector, nbands  )
            E2=wannier_bands(HWannier, cell_map, kvector+qnormalized, nbands  )
            
            V1=wannier_vectors(HWannier, cell_map, kvector)
            V2=wannier_vectors(HWannier, cell_map, kvector+qnormalized )

            for lower in 1:valence_bands+1
                for upper in valence_bands+1:nbands
                    Elower = E1[lower]
                    Eupper = E2[upper]
                    overlap=(np.abs(np.dot(V1[:, lower], np.conj(V2[:, upper]))))^2;
                    f1=np.heaviside( μ-Elower, 0.5)
                    f2=np.heaviside( μ-Eupper, 0.5)

                    DeltaE=Eupper-Elower
                    if DeltaE>0
                        Polarization_Array[round(Int, histogram_width*DeltaE+1)] = Polarization_Array[round(Int, histogram_width*DeltaE+1)]+π*(f2-f1)/V*overlap*(1/mesh)^2*histogram_width*spin
                    end
                end
            end
        end
    end

    return Polarization_Array
end


function im_polarization(wannier_file_up::String, wannier_file_dn::String,  cell_map_file_up::String, cell_map_file_dn::String, nbands::Int, valence_bands_up::Int, valence_bands_dn::Int, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; kwargs...) 
    #Here we add the independent polarizations from different spin channels 
    spin_up_pol = im_polarization(wannier_file_up, cell_map_file_up, nbands, valence_bands_up, lattice_vectors, q, μ; kwargs... )
    spin_dn_pol = im_polarization(wannier_file_dn, cell_map_file_dn, nbands, valence_bands_dn, lattice_vectors, q, μ; kwargs... )
    return (spin_up_pol + spin_dn_pol)
end

function im_polarization(HWannierup::Array{Float64, 3}, HWannierdn::Array{Float64, 3},  cell_map_up::Array{Float64, 2}, cell_map_dn::Array{Float64, 2}, nbands::Int, valence_bands_up::Int, valence_bands_dn::Int, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, μ::Real; kwargs...)
    #Here we add the independent polarizations from different spin channels 
    spin_up_pol = im_polarization(HWannierup, cell_map_up, nbands, valence_bands_up, lattice_vectors, q, μ; kwargs... )
    spin_dn_pol = im_polarization(HWannierdn, cell_map_dn, nbands, valence_bands_dn, lattice_vectors, q, μ; kwargs... )
    return (spin_up_pol + spin_dn_pol)
end

"Applies the kramers-kronig relations onto a 1 dimensional array of numbers consisting of the imaginary value of the polarization to return the real value of polarization"
function kramers_kronig(ω::Real, im_pol::Array{<:Real, 1}, max_energy::Real, histogram_width::Real) 
    omegaprime=collect(1:histogram_width*max_energy).*1/histogram_width
    sum(1/histogram_width*2/π*im_pol.*omegaprime./(omegaprime.^2 .- (ω+ω*0.03im)^2))

end

"Mostly provided for checking the reciprocity relation between the real and imaginary susceptibilities. Give the real susceptibility to obtain the imaginary susceptibility"
function kramers_kronig_reverse(ω::Real, re_pol::Array{<:Real, 1}, max_energy::Real, domega::Real) 

    omegaprime=collect(0:domega:max_energy)
    sum(-domega*2/π*re_pol.*ω./(omegaprime.^2 .- (ω+ω*0.03im)^2))

end

function kramers_kronig_reverse_scipy(ω::Real, re_pol::Array{<:Real, 1}, max_energy::Real, domega::Real, max_energy_integration::Real; kwargs...) 
    interpolated_res=interpol.interp1d(0:domega:max_energy, re_pol)
    
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=-2/pi*interpolated_res(omegaprime)*ω/(omegaprime+ω)

    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]

end

"Applies the kramers-kronig relations but with scipy's cauchy weight; kwargs for scipy.integrate.quad supported"
function kramers_kronig_scipy(ω::Real, im_pol::Array{<:Real, 1}, max_energy::Real, histogram_width::Real, max_energy_integration::Real; kwargs...) 
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(max_energy-1/histogram_width), im_pol)
    
    ErrorAbs=1e-20
    cauchy_inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime+ω)

    return pyintegrate.quad(cauchy_inner_function, 0, max_energy_integration, weight="cauchy",  epsrel=ErrorAbs, epsabs=ErrorAbs, limit=75,  wvar= ω ; kwargs...)[1]

end

function kramers_kronig_reverse_quadgk(ω::Real, re_pol::Array{<:Real, 1}, max_energy::Real, domega::Real, max_energy_integration::Real ; δ::Real = 0.1, kwargs...) 
    
    interpolated_res=interpol.interp1d(0:domega:max_energy, re_pol)
    
    inner_function(omegaprime)=-2/pi*interpolated_res(omegaprime)*ω/(omegaprime^2-(ω+1im*δ)^2)

    return real(quadgk(inner_function, 0, max_energy_integration; kwargs...)[1])

end

function kramers_kronig_quadgk(ω::Real, im_pol::Array{<:Real, 1}, max_energy::Real, histogram_width::Real, max_energy_integration::Real; δ::Real = 0.1, kwargs...) 
    
    interpolated_ims=interpol.interp1d(0:1/histogram_width:(max_energy-1/histogram_width), im_pol)
    
    inner_function(omegaprime)=2/pi*interpolated_ims(omegaprime)*omegaprime/(omegaprime^2-(ω+1im*δ)^2)

    return real(quadgk(inner_function, 0, max_energy_integration; kwargs...)[1])

end

function epsilon_integrand(wannier_file::String, cell_map_file::String, k₁::Real, k₂::Real, q::Array{<:Real, 1}, μ::Real, ω::Real, ϵ::Real; spin::Int=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(wannier_file, cell_map_file, kvector,  )
    ϵ₂ =wannier_bands(wannier_file, cell_map_file, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    real(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

function epsilon_integrand_imaginary(wannier_file::String, cell_map_file::String, k₁::Real, k₂::Real, q::Array{<:Real, 1}, μ::Real, ω::Real, ϵ::Real; spin::Int=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(wannier_file, cell_map_file, kvector  )
    ϵ₂ =wannier_bands(wannier_file, cell_map_file, kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    imag(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

function epsilon_integrand(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k₁::Real, k₂::Real, q::Array{<:Real, 1}, μ::Real, ω::Real, ϵ::Real; spin::Int=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(HWannier, cell_map, kvector  )
    ϵ₂ =wannier_bands(HWannier, cell_map,  kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    real(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end

function epsilon_integrand_imaginary(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, k₁::Real, k₂::Real, q::Array{<:Real, 1}, μ::Real, ω::Real, ϵ::Real; spin::Int=1)
    kvector=[k₁, k₂, 0]
    ϵ₁ =wannier_bands(HWannier, cell_map,  kvector  )
    ϵ₂ =wannier_bands(HWannier, cell_map,  kvector+q  )
    f = ϵ₁<μ ? 1 : 0
    imag(1/(2π)^2*spin*2*f*(ϵ₁-ϵ₂)/((ϵ₁-ϵ₂)^2-(ω+1im*ϵ)^2))
end


function direct_epsilon(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, ω::Real, μ::Real; spin::Int = 1, ϵ::Real = 0.01, kwargs...) 
    
    kwargsdict=Dict()

    for kwarg in kwargs
        push!(kwargsdict, kwarg.first => kwarg.second)
    end

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*pyintegrate.nquad((k₁, k₂) -> epsilon_integrand(wannier_file, cell_map_file, k₁, k₂, qnormalized, μ, ω, ϵ, spin=spin), [[0, 1], [0, 1]], opts=kwargsdict)[1]

    1-e²ϵ/(2qabs)*polarization

end

function direct_epsilon(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, ω::Real, μ::Real; spin::Int = 1, ϵ::Real = 0.01, kwargs...) 
    
    kwargsdict=Dict()

    for kwarg in kwargs
        push!(kwargsdict, kwarg.first => kwarg.second)
    end

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*pyintegrate.nquad((k₁, k₂) -> epsilon_integrand(HWannier, cell_map, k₁, k₂, qnormalized, μ, ω, ϵ, spin=spin), [[0, 1], [0, 1]], opts=kwargsdict)[1]

    1-e²ϵ/(2qabs)*polarization

end


"Direct 2D integration for Epsilon with HCubature"
function direct_epsilon_cubature(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, ω::Real, μ::Real; spin::Int = 1, ϵ::Real = 0.01, kwargs...)

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand(wannier_file, cell_map_file, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]

    1-e²ϵ/(2qabs)*polarization

end

function direct_epsilon_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, ω::Real, μ::Real; spin::Int = 1, ϵ::Real = 0.01, kwargs...)

    qnormalized = normalize_kvector(lattice_vectors, q)
    qabs=sqrt(sum(q.^2))

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand(HWannier, cell_map, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]

    1-e²ϵ/(2qabs)*polarization

end


"Find the imaginary value of polarization through hcubature "
function im_polarization_cubature(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, ω::Real, μ::Real; spin::Int=1, ϵ::Real=0.01, kwargs...) 

    qnormalized = normalize_kvector(lattice_vectors, q)

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand_imaginary(wannier_file, cell_map_file, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]

    return polarization

end

function im_polarization_cubature(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, q::Array{<:Real, 1}, ω::Real, μ::Real; spin::Int=1, ϵ::Real=0.01, kwargs...) 

    qnormalized = normalize_kvector(lattice_vectors, q)

    brillouin_area=brillouin_zone_area(lattice_vectors) 
    
    polarization=brillouin_area*hcubature((k) -> epsilon_integrand_imaginary(HWannier, cell_map, k[1], k[2], qnormalized, μ, ω, ϵ, spin=spin), [0, 0], [1, 1]; kwargs...)[1]

    return polarization

end


"returns the non-local, non-static dielectric function"
function return_2d_epsilon(q::Real, ω::Real, im_pol::Array{<:Real, 1}, max_energy::Real, histogram_width::Real) 
    return 1-e²ϵ/abs(2q)*kramers_kronig(ω, im_pol, max_energy, histogram_width)
end

"returns the non-local, non-static dielectric function using scipy functionality"
function return_2d_epsilon_scipy(q::Real, ω::Real, im_pol::Array{<:Real, 1}, max_energy::Real, histogram_width::Real, max_energy_integration::Real) 
    return 1-e²ϵ/abs(2q)*kramers_kronig_scipy(ω, im_pol, max_energy, histogram_width, max_energy_integration)
end

function return_2d_epsilon_quadgk(q::Real, ω::Real, im_pol::Array{<:Real, 1}, max_energy::Real, histogram_width::Real, max_energy_integration::Real; δ::Real = 0.1, kwargs... )
    return 1-e²ϵ/(2abs(q))*kramers_kronig_quadgk(ω, im_pol, max_energy, histogram_width, max_energy_integration; δ, kwargs...)  
end