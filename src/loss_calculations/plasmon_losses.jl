
#We calculate losses at 0th (landau damping), 1st, and 2nd orders in phonon-assisted damping

function landau_damping(wannier_file::String, cell_map_file::String, histogram_length::Int, mesh::Int, q::Array{T, 1}, μ::R, offset, energy_range) where {T<:Number, R<:Number}
    lossarray = zeros(historam_width*energy_range)
    qabs = sqrt(sum(q.^2))
    qnormalized = normalize_kvector(q)
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

function first_order_damping()
    return 2π/ħ
end

function second_order_damping()
    return 2π/ħ
end

