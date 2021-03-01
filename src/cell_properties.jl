using Plots
using PyCall
using LinearAlgebra 
using Distances

function cell_vectors(lattice_file::String)
    run(`cat $lattice_file`);
    run(`pwd`)
end

"reciprocal_vectors returns the reciprocal lattice vectors when supplied with three real space vectors"
function reciprocal_vectors(lattice_vectors::Array{Array{T, 1},1}) where T <: Number
    
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]

    V=dot(a1, cross(a2, a3))
    b1=2π/V*cross(a2, a3)
    b2=2π/V*cross(a3, a1)
    b3=2π/V*cross(a1, a2)
    return b1, b2, b3

end

function in_wigner_seitz(lattice_vectors::Array{Array{T, 1},1}, rvec::Array{Array{R, 1}, 1}) where {T<:Number, R<:Number}
    
    vec1 = lattice_vectors[1]
    vec2 = lattice_vectors[2]
    vec3 = lattice_vectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            for k in -2:-2
                current_vec = vec1*i+vec2*j+vec3*k
                push!(distances_array, euclidean(currentvec, rvec) )
            end
        end
    end

    if euclidean(rvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end

function in_wigner_seitz(lattice_vectors::lattice, rvec::Array{Array{R, 1}, 1}) where {T<:Number, R<:Number}
    
    vec1 = lattice_vectors.rvectors[:, 1]*bohrtoangstromn
    vec2 = lattice_vectors.rvectors[:, 2]*bohrtoangstromn
    vec3 = lattice_vectors.rvectors[:, 3]*bohrtoangstromn
    distances_array = []
    for i in -2:2
        for j in -2:2
            for k in -2:-2
                current_vec = vec1*i+vec2*j+vec3*k
                push!(distances_array, euclidean(currentvec, rvec) )
            end
        end
    end

    if euclidean(rvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end


function in_brillouin(lattice_vectors::Array{Array{T, 1},1}, kvec::Array{Array{R, 1}, 1}) where {T<:Number, R<:Number}
    
    bvectors = reciprocal_vectors(lattice_vectors)

    vec1 = bvectors[1]
    vec2 = bvectors[2]
    vec3 = bvectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            for k in -2:-2
                current_vec = vec1*i+vec2*j+vec3*k
                push!(distances_array, euclidean(currentvec, kvec) )
            end
        end
    end

    if euclidean(kvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end

function in_brillouin(lattice_vectors::lattice, kvec::Array{Array{R, 1}, 1}) where {T<:Number, R<:Number}
    
    lattice_vectors_array = [lattice_vectors.rvectors[:, 1]*bohrtoangstromn,lattice_vectors.rvectors[:, 2]*bohrtoangstromn, lattice_vectors.rvectors[:, 3]*bohrtoangstromn ]
    bvectors = reciprocal_vectors(lattice_vectors_array)

    vec1 = bvectors[1]
    vec2 = bvectors[2]
    vec3 = bvectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            for k in -2:-2
                current_vec = vec1*i+vec2*j+vec3*k
                push!(distances_array, euclidean(currentvec, kvec) )
            end
        end
    end

    if euclidean(kvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end



"returns the normalized kvector (in the basis of the reciprocal lattice vectors)"
function normalize_kvector(lattice_vectors::Array{Array{T, 1},1}, unnormalized_kvector) where T <: Number

    b1, b2, b3 = reciprocal_vectors(lattice_vectors)

    vectors_array=Array{Float64,2}(undef, (3, 3))
    
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3

    inv(vectors_array)*unnormalized_kvector

end

function unnormalize_kvector(lattice_vectors::Array{Array{T, 1},1}, normalized_kvector) where T <: Number

    b1, b2, b3 = reciprocal_vectors(lattice_vectors)

    vectors_array=Array{Float64,2}(undef, (3, 3))
    
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3

    vectors_array*normalized_kvector

end


"Returns the 2d area of the lattice. The assumption is made that the lattice is in the x-y plane"
function brillouin_zone_area(lattice_vectors::Array{Array{T, 1},1}) where T <: Number

    b_vectors=reciprocal_vectors(lattice_vectors)

    b_vectors_2d = []
    for b_vector in b_vectors 
        if b_vector[3] ≈ 0
            push!(b_vectors_2d, b_vector)

        end
    end

    b2d_1, b2d_2 = b_vectors_2d 

    return sqrt(sum(cross(b2d_1, b2d_2).^2))

end

function ion_positions(ionpos_file::String)
    run(`cat $ionpos_file`);
    run(`pwd`)
end

function plot_lattice(lattice_file::String)
    ion_position_vectors=String[]
    open(lattice_file, "r") do io
        readline(io)
        ion_position_vectors=readlines(io);
    end
end

function pyversion()
    sys=pyimport("sys")
    print("You are currently running this version of python: $(sys.executable)")
end
