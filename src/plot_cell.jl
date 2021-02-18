using Plots
using PyCall
using LinearAlgebra 
using Distances

function cell_vectors(lattice_file::String)
    run(`cat $lattice_file`);
    run(`pwd`)
end

function density_of_states(dosfile_1::String, dosfile_2::String )
   
    np=pyimport("numpy")
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up")
    plot!(np.loadtxt(dosfile_2)[:, 1]*27.2, np.loadtxt(dosfile_2)[:, 2]/27.2, linewidth=4,  size=(800, 400), label="Spin Down")

end

function density_of_states(dosfile_1::String)

    np=pyimport("numpy")
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized")

end

function density_of_states_wannier(wannier_file::String, cell_map_file::String; mesh=100, histogram_width=100, energy_range=10, offset=0)
    np=pyimport("numpy")
    WannierDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ=  wannier_bands(wannier_file, cell_map_file, [x_mesh/mesh, y_mesh/mesh, 0])
            WannierDOS[round(Int, histogram_width*(ϵ+offset))]=WannierDOS[round(Int, histogram_width*(ϵ+offset))]+histogram_width*(1/mesh)^2

        end
    end

    return WannierDOS

end

function density_of_states_wannier(wannier_file_up::String, cell_map_file_up::String, wannier_file_dn::String, cell_map_file_dn::String,; mesh=100, histogram_width=100, energy_range=10, offset=0)
   
    np=pyimport("numpy")
    WannierDOSUp=np.zeros(histogram_width*energy_range)
    WannierDOSDn=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ_up=  wannier_bands(wannier_file_up, cell_map_file_up, [x_mesh/mesh, y_mesh/mesh, 0])
            WannierDOSUp[round(Int, histogram_width*(ϵ_up+offset))]=WannierDOSUp[round(Int, histogram_width*(ϵ_up+offset))]+histogram_width*(1/mesh)^2

            ϵ_dn=  wannier_bands(wannier_file_dn, cell_map_file_dn, [x_mesh/mesh, y_mesh/mesh, 0])
            WannierDOSDn[round(Int, histogram_width*(ϵ_dn+offset))]=WannierDOSDn[round(Int, histogram_width*(ϵ_dn+offset))]+histogram_width*(1/mesh)^2

        end
    end

    return WannierDOSUp, WannierDOSDn

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

function plot_bands(band_file::String, num_bands::Int, num_points::Int)

    reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
    exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*27.2;
    exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*27.2;

    plot(exactenergiesdown, color="black", label="", linewidth=2)
    plot!(exactenergiesup, color="purple", label="", linewidth=2)

end

function plot_phonons(cell_map::String, phononOmegaSq::String, kpoints::String)
    np=pyimport("numpy")
    cellMapPh = np.loadtxt(cell_map)[:,1:3]
    forceMatrixPh = np.fromfile(phononOmegaSq, dtype=np.float64)
    nCellsPh = size(cellMapPh)[1]
    nModesPh = Int(np.sqrt(size(forceMatrixPh)[1] / nCellsPh))
    forceMatrixPh = np.reshape(forceMatrixPh, (nCellsPh,nModesPh,nModesPh))

    eV=1/27.2
    kpointsIn = np.loadtxt(kpoints, skiprows=2, usecols=(1,2,3))
    nKin = size(kpointsIn)[1]
    #--- Fourier transform from real to k space:
    forceMatrixTilde = np.tensordot(np.exp((2im*np.pi)*np.dot(kpointsIn,transpose(cellMapPh))), forceMatrixPh, axes=1)
    #--- Diagonalize:
    omegaSq, normalModes = np.linalg.eigh(forceMatrixTilde)
    plot(title="Phonon Dispersion", titlefontsize=20, ytickfontsize=15, sqrt.(abs.(omegaSq))/eV, linewidth=2, color="orange", legend=false, size=(800, 1000), xticks=[])

end

function pyversion()
    sys=pyimport("sys")
    print("You are currently running this version of python: $(sys.executable)")
end

function kramers_kronig(ω::T, im_pol::Array{R, 1}, max_energy::S, histogram_width::Q) where {T<:Number, R<:Number, Q<:Number, S<:Number}
    #pyintegrate=pyimport("scipy.integrate")

    #pyintegrate.quad(sin, 0, max_energy)
    #ω_array=collect(1/histogram_width:1/histogram_width:max_energy)
    #sum(im_pol.*1 ./ ( ω_array.-ω.+0.001*im))*1/histogram_width
    omegaprime=collect(1:histogram_width*max_energy).*1/histogram_width
    sum(1/histogram_width*2/π*im_pol.*omegaprime./(omegaprime.^2 .- (ω+ω*0.03im)^2))

end


function wannier_bands(wannier_file::String, cell_map_file::String, k::Array{T, 1}) where T<:Number
    np=pyimport("numpy")
    cell_map=np.loadtxt(cell_map_file)
    Hwannier=permutedims(reshape(np.loadtxt(wannier_file), (43, 1, 1)), [1, 3, 2])
    eV = 1/27.2114 #in Hartrees
    phase = np.exp(2im*np.pi*cell_map*k); H = np.tensordot(phase, Hwannier, axes=1); E, U=np.linalg.eigh(H);
    return E[1]/eV 
end