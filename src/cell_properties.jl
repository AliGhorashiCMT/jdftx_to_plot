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

function in_wigner_seitz(a1::Array{T, 1}, a2::Array{T, 1}) where T<:Number
    return euclidean(a1, a2)
end

function normalize_kvector(lattice_vectors::Array{Array{T, 1},1}, unnormalized_kvector) where T <: Number

    b1, b2, b3 = reciprocal_vectors(lattice_vectors)

    vectors_array=Array{Float64,2}(undef, (3, 3))
    
    vectors_array[:, 1], vectors_array[:, 2], vectors_array[:, 3] = b1, b2, b3

    inv(vectors_array)*unnormalized_kvector

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
