function cell_vectors(lattice_file::String)
    run(`cat $lattice_file`);
    run(`pwd`)
end

"reciprocal_vectors returns the reciprocal lattice vectors when supplied with three real space vectors"
function reciprocal_vectors(lattice_vectors::Array{<:Array{<:Real, 1},1}) 
    
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]

    V=dot(a1, cross(a2, a3))
    b1=2π/V*cross(a2, a3)
    b2=2π/V*cross(a3, a1)
    b3=2π/V*cross(a1, a2)
    return b1, b2, b3

end

function in_wigner_seitz(lattice_vectors::Array{<:Array{<:Real, 1},1}, rvec::Array{<:Real, 1}) 
    vec1 = lattice_vectors[1]
    vec2 = lattice_vectors[2]
    vec3 = lattice_vectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, rvec) )
            end
        end
    end

    if euclidean(rvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end

function in_wigner_seitz(lattice_vectors::lattice, rvec::Array{<:Real, 1}) 
    
    vec1 = lattice_vectors.rvectors[:, 1]*bohrtoangstrom
    vec2 = lattice_vectors.rvectors[:, 2]*bohrtoangstrom
    vec3 = lattice_vectors.rvectors[:, 3]*bohrtoangstrom
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
                
            else 
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, rvec) )
            end
        end
    end

    if euclidean(rvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end



function in_brillouin(lattice_vectors::Array{<:Array{<:Real, 1},1}, kvec::Array{<:Real, 1}) 
    
    bvectors = reciprocal_vectors(lattice_vectors)

    vec1 = bvectors[1]
    vec2 = bvectors[2]
    vec3 = bvectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else 
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, kvec) )
            end
        end
    end

    if euclidean(kvec, [0, 0, 0]) < minimum(distances_array)
        return true
    else 
        return false
    end

end

function in_brillouin(lattice_vectors::lattice, kvec::Array{<:Real, 1}) 
    
    lattice_vectors_array = [lattice_vectors.rvectors[:, 1]*bohrtoangstrom,lattice_vectors.rvectors[:, 2]*bohrtoangstrom, lattice_vectors.rvectors[:, 3]*bohrtoangstrom ]
    bvectors = reciprocal_vectors(lattice_vectors_array)

    vec1 = bvectors[1]
    vec2 = bvectors[2]
    vec3 = bvectors[3]
    distances_array = []
    for i in -2:2
        for j in -2:2
            if i==0 && j==0
                continue
            else 
                current_vec = vec1*i+vec2*j
                push!(distances_array, euclidean(current_vec, kvec) )
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
function normalize_kvector(lattice_vectors::Array{<:Array{<:Real, 1},1}, unnormalized_kvector)

    b1, b2, b3 = reciprocal_vectors(lattice_vectors)

    vectors_array=Array{Float64,2}(undef, (3, 3))
    
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3

    inv(vectors_array)*unnormalized_kvector

end

function unnormalize_kvector(lattice_vectors::Array{<:Array{<:Real, 1},1}, normalized_kvector) 

    b1, b2, b3 = reciprocal_vectors(lattice_vectors)

    vectors_array=Array{Float64,2}(undef, (3, 3))
    
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3

    vectors_array*normalized_kvector

end

function unit_cell_area(lattice_vectors::Array{<:Array{<:Real, 1},1}) 
    a1, a2, a3 = lattice_vectors[1], lattice_vectors[2], lattice_vectors[3]

    A=sqrt(dot(cross(a1, a2), cross(a1, a2)))
    return A

end

function unit_cell_area(lattice_vectors::lattice)
    a1, a2, a3 = lattice_vectors.rvectors[:, 1]*bohrtoangstrom, lattice_vectors.rvectors[:, 2]*bohrtoangstrom, lattice_vectors.rvectors[:, 3]*bohrtoangstrom

    A=sqrt(dot(cross(a1, a2), cross(a1, a2)))
    return A

end


"Returns the 2d area of the lattice. The assumption is made that the lattice is in the x-y plane"
function brillouin_zone_area(lattice_vectors::Array{<:Array{<:Real, 1},1}) 

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
