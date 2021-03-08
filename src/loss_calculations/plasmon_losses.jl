#We calculate losses at 0th (landau damping), 1st, and 2nd orders in phonon-assisted damping

function landau_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1},1}, histogram_length::Int, mesh::Int, q::Array{<:Real, 1}, μ::Real, energy_range::Real) 
    lossarray = zeros(histogram_width*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵ1 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0])
            ϵ2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
            f1 = ϵ1<μ ? 1 : 0
            f2 = ϵ2>μ ? 1 : 0
            if f1>0 && f2>0
                ω = ϵ2-ϵ1
                lossarray[round(Int, (ω+offset)*histogram_length  )] = lossarray[round(Int, (ω+offset)*histogram_length  )] + 2π/ħ*e²ϵ/4*ω/qabs*f1*f2*(1/mesh)^2*histogram_length
            end
        end
    end
    return lossarray
end

function landau_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1},1}, histogram_length::Int, mesh::Int, q::Array{<:Real, 1}, μ::Real, energy_range::Real) 
    lossarray = zeros(histogram_width*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵ1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            ϵ2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)
            f1 = ϵ1<μ ? 1 : 0
            f2 = ϵ2>μ ? 1 : 0
            if f1>0 && f2>0
                ω = ϵ2-ϵ1
                lossarray[round(Int, (ω+offset)*histogram_length  )] = lossarray[round(Int, (ω+offset)*histogram_length  )] + 2π/ħ*e²ϵ/4*ω/qabs*f1*f2*(1/mesh)^2*histogram_length
            end
        end
    end
    return lossarray
end

function first_order_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)
    for xmesh in 1:mesh
        for ymesh in 1:mesh
            ϵinitial = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0])
            ϵmiddle = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0
            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    ϵmiddle2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    ϵfinal = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                    ffinal = ϵfinal>μ ? 1 : 0

                    ω = ϵfinal-ϵinitial+ϵphonon
                    if ω>0
                        lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*(1/(ϵmiddle-ϵinitial-ω)+1/(ϵmiddle2-ϵinitial+ϵphonon))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*fmiddle1*ffinal*(1/mesh)^4*histogram_length
                    end
                end
            end
        end
    end
    return lossarray
end



function second_order_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    

    #= In all cases, we consider an initial state with an electron at kx, ky = xmesh/mesh, ymesh/mesh and a plasmon with wavevector q 
     Furthermore, the final state is an electron with momentum (xmesh+xmesh1+xmesh2)/mesh, (ymesh+ymesh1+ymesh2)/mesh - qnormalized
    as well as two emitted phonons. We sum over all intermmediate states and square the sum in the integrand. The three possible decay channels are 
    plasmon absorption first, second, or third. Note that thus for each iteration of the sum, we're summing contributions from virtual 
    electronic intermmediate states. Note also that the matrix elements included here contain Fermi occupation functions 
     =#

    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)

    for xmesh in 1:mesh
        for ymesh in 1:mesh

            ϵinitial = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0])
            #first middle state (plasmon absorbed first)
            ϵmiddle = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0

            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    #second middle state (phonon absorbed first)
                    ϵmiddle2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    #Then consider second middle states from one phonon and plasmon absoption first, (phonon absorbed last)

                    ϵsecondmiddle2 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0] + qnormalized)        
                    fsecondmiddle2 = ϵsecondmiddle1>μ ? 1 : 0


                    for xmesh2 in 1:mesh
                        for ymesh2 in 1:mesh

                            #first consider second middle states originating from 2 phonons being absorbed first (plasmon absorbed last)
                            ϵsecondmiddle1 = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0])        
                            fsecondmiddle1 = ϵsecondmiddle1>μ ? 1 : 0

                            #The final energy will always be that of the electronic state corresponding to the original kvector plus the two phonon kvectors and the plasmon kvector
                            ϵfinal = wannier_bands(wannier_file, cell_map_file, [xmesh/mesh, ymesh/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                            ffinal = ϵfinal>μ ? 1 : 0

                            ω = ϵfinal-ϵinitial+ϵphonon
                            if ω>0
                                lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*( 1/(ϵmiddle-ϵinitial-ω)*1/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon)*fmiddle1*fsecondmiddle2 + fmiddle2/(ϵmiddle2-ϵinitial+ϵphonon)*( fsecondmiddle1/(ϵsecondmiddle1-ϵinitial+2*ϵphonon) + fsecondmiddle2/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon) ))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^6*histogram_length
                            end
                        end
                    end
                end
            end
        end
    end
    return lossarray
end


function second_order_damping(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}, q::Array{<:Real, 1}, μ::Real, ϵphonon::Real, gph::Real; histogram_length::Real=100, mesh::Int=30, energy_range::Real=10) 
    

    #= In all cases, we consider an initial state with an electron at kx, ky = xmesh/mesh, ymesh/mesh and a plasmon with wavevector q 
     Furthermore, the final state is an electron with momentum (xmesh+xmesh1+xmesh2)/mesh, (ymesh+ymesh1+ymesh2)/mesh - qnormalized
    as well as two emitted phonons. We sum over all intermmediate states and square the sum in the integrand. The three possible decay channels are 
    plasmon absorption first, second, or third. Note that thus for each iteration of the sum, we're summing contributions from virtual 
    electronic intermmediate states. Note also that the matrix elements included here contain Fermi occupation functions 
     =#

    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
    cell_area = unit_cell_area(lattice_vectors)

    for xmesh in 1:mesh
        for ymesh in 1:mesh

            ϵinitial = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0])
            #first middle state (plasmon absorbed first)
            ϵmiddle = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+qnormalized)        

            finitial = ϵinitial<μ ? 1 : 0
            fmiddle1 = ϵmiddle>μ ? 1 : 0

            for xmesh1 in 1:mesh
                for ymesh1 in 1:mesh

                    #second middle state (phonon absorbed first)
                    ϵmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0])        
                    fmiddle2 = ϵmiddle2>μ ? 1 : 0

                    #Then consider second middle states from one phonon and plasmon absoption first, (phonon absorbed last)

                    ϵsecondmiddle2 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0] + qnormalized)        
                    fsecondmiddle2 = ϵsecondmiddle1>μ ? 1 : 0


                    for xmesh2 in 1:mesh
                        for ymesh2 in 1:mesh

                            #first consider second middle states originating from 2 phonons being absorbed first (plasmon absorbed last)
                            ϵsecondmiddle1 = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh1/mesh, ymesh1/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0])        
                            fsecondmiddle1 = ϵsecondmiddle1>μ ? 1 : 0

                            #The final energy will always be that of the electronic state corresponding to the original kvector plus the two phonon kvectors and the plasmon kvector
                            ϵfinal = wannier_bands(HWannier, cell_map, [xmesh/mesh, ymesh/mesh, 0]+[xmesh2/mesh, ymesh2/mesh, 0]+qnormalized+[xmesh1/mesh, ymesh1/mesh, 0])        
                    
                            ffinal = ϵfinal>μ ? 1 : 0

                            ω = ϵfinal-ϵinitial+ϵphonon
                            if ω>0
                                lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + 1/cell_area*( 1/(ϵmiddle-ϵinitial-ω)*1/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon)*fmiddle1*fsecondmiddle2 + fmiddle2/(ϵmiddle2-ϵinitial+ϵphonon)*( fsecondmiddle1/(ϵsecondmiddle1-ϵinitial+2*ϵphonon) + fsecondmiddle2/(ϵsecondmiddle2-ϵinitial-ω+ϵphonon) ))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*ffinal*(1/mesh)^6*histogram_length
                            end
                        end
                    end
                end
            end
        end
    end
    return lossarray
end



