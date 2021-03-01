
#We calculate losses at 0th (landau damping), 1st, and 2nd orders in phonon-assisted damping

function landau_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{Array{S, 1},1}, histogram_length::Int, mesh::Int, q::Array{T, 1}, μ::R, offset, energy_range) where {T<:Number, R<:Number, S<:Number}
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

function first_order_damping(wannier_file::String, cell_map_file::String, lattice_vectors::Array{Array{S, 1}, 1}, q::Array{T, 1}, μ::R, ϵphonon, gph; histogram_length=100, mesh=30, energy_range=10) where {T<:Number, R<:Number, S<:Number}
    lossarray = zeros(histogram_length*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(lattice_vectors, q)
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
                        lossarray[round(Int, ω*histogram_length+1)] = lossarray[round(Int, ω*histogram_length + 1 )] + (1/(ϵmiddle-ϵinitial-ω)+1/(ϵmiddle2-ϵinitial+ϵphonon))^2*gph^2*2π/ħ*e²ϵ/4*ω/qabs*finitial*fmiddle1*ffinal*(1/mesh)^4*histogram_length
                    end
                end
            end
        end
    end
    return lossarray
end



function second_order_damping()
    return 2π/ħ
end

