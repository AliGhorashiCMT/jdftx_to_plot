function density_of_states(dosfile_1::String, dosfile_2::String; kwargs... )
   
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    plot!(np.loadtxt(dosfile_2)[:, 1]*27.2, np.loadtxt(dosfile_2)[:, 2]/27.2, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)

end

function density_of_states(dosfile_1::String; kwargs...)

    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)

end

"""
The typical density of states outputed by JDFTX is per unit cell. However, sometimes it is more relevant to know the 
density of states per unit volume. This is simply equivalent to dividing the conventional DOS by the unit cell 
volume 
"""
function density_of_states_per_area(dosfile_1::String, lattice_vecs::Array{<:Array{<:Real, 1}, 1}; kwargs...)

    ucell_area = unit_cell_area(lattice_vecs)
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)

end

function density_of_states_per_area(dosfile_1::String, lattice_vecs::lattice; kwargs...)

    ucell_area = unit_cell_area(lattice_vecs)
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)

end


function density_of_states_per_area(dosfile_1::String, dosfile_2::String, lattice_vecs::Array{<:Array{<:Real, 1}, 1}; kwargs... )
    ucell_area = unit_cell_area(lattice_vecs)
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    plot!(np.loadtxt(dosfile_2)[:, 1]*27.2, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]/27.2, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)

end

function density_of_states_per_area(dosfile_1::String, dosfile_2::String, lattice_vecs::lattice; kwargs... )
    ucell_area = unit_cell_area(lattice_vecs)
    plot(np.loadtxt(dosfile_1)[:, 1]*27.2, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    plot!(np.loadtxt(dosfile_2)[:, 1]*27.2, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]/27.2, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)

end

function density_of_states_wannier_quad(wannier_file::String, cell_map_file::String, ϵ::Real; δ=.1, kwargs...) 

    1/π*hcubature(vec->imag(-1/(ϵ-wannier_bands(wannier_file, cell_map_file, [vec[1], vec[2], 0])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]

end

function density_of_states_wannier_quad(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵ::Real; δ::Real = 0.1, kwargs...)

    1/π*hcubature(vec->imag(-1/(ϵ-wannier_bands(HWannier, cell_map, [vec[1], vec[2], 0])+1im*δ)), [0, 0], [1, 1]; kwargs...)[1]

end

function density_of_states_wannier_scipy_quad(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵ::Real; δ::Real = 0.1, kwargs...) 

    nquad = pyintegrate.nquad

    optdict=Dict()

    for kwarg in kwargs
        push!(optdict, kwarg[1]=>kwarg[2])
    end

    1/π*nquad((x, y)->imag(-1/(ϵ-wannier_bands(HWannier, cell_map, [x, y, 0])+1im*δ)), [[0, 1], [0, 1]], opts=optdict)[1]

end



function density_of_states_wannier_scipy_quad(wannier_file::String, cell_map_file::String, ϵ::Real; δ::Real = 0.1, kwargs...) 

    nquad = pyintegrate.nquad

    optdict=Dict()

    for kwarg in kwargs
        push!(optdict, kwarg[1]=>kwarg[2])
    end

    1/π*nquad((x, y)->imag(-1/(ϵ-wannier_bands(wannier_file, cell_map_file, [x, y, 0])+1im*δ)), [[0, 1], [0, 1]], opts=optdict)[1]

end


function density_of_states_wannier_quad_check(wannier_file::String, cell_map_file::String, ϵmin::Real, ϵmax::Real, numpoints::Int; δ=.1, kwargs...) 

    ϵdif=(ϵmax-ϵmin)/numpoints
    dosarray=[]
    for i in 0:numpoints
        ϵ=ϵmin+ϵdif*i
        push!(dosarray, density_of_states_wannier_quad(wannier_file, cell_map_file, ϵ; δ, kwargs... ))
    end
    return sum(ϵdif*dosarray)

end

function density_of_states_wannier_quad_check(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, ϵmin::Real, ϵmax::Real, numpoints::Int; δ::Real = 0.1, kwargs...) 

    ϵdif=(ϵmax-ϵmin)/numpoints
    dosarray=[]
    for i in 0:numpoints
        ϵ=ϵmin+ϵdif*i
        push!(dosarray, density_of_states_wannier_quad(HWannier, cell_map, ϵ; δ, kwargs... ))
    end
    return sum(ϵdif*dosarray)

end


function density_of_states_wannier(wannier_file::String, cell_map_file::String; mesh::Int = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ=  wannier_bands(wannier_file, cell_map_file, [x_mesh/mesh, y_mesh/mesh, 0])
            WannierDOS[round(Int, histogram_width*(ϵ+offset))]=WannierDOS[round(Int, histogram_width*(ϵ+offset))]+histogram_width*(1/mesh)^2

        end
    end

    return WannierDOS

end

function density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}; mesh::Int = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ=  wannier_bands(HWannier, cell_map, [x_mesh/mesh, y_mesh/mesh, 0])
            WannierDOS[round(Int, histogram_width*(ϵ+offset))]=WannierDOS[round(Int, histogram_width*(ϵ+offset))]+histogram_width*(1/mesh)^2

        end
    end

    return WannierDOS

end

"""
The standard DOS function using wannier functions returns the density of states per eV per unit cell. 
At times it is more convenient to obtain the DOS per eV per angstrom^2 
"""
function density_of_states_wannier_per_area(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, lattice_vectors::Array{<:Array{<:Real, 1}, 1}; mesh::Int = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    ucell_area = unit_cell_area(lattice_vectors)
    WannierDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ=  wannier_bands(HWannier, cell_map, [x_mesh/mesh, y_mesh/mesh, 0])
            WannierDOS[round(Int, histogram_width*(ϵ+offset))]=WannierDOS[round(Int, histogram_width*(ϵ+offset))]+histogram_width*(1/mesh)^2

        end
    end

    return WannierDOS/ucell_area

end


function density_of_states_wannier(wannier_file::String, cell_map_file::String, nbands::Int; exclude_bands::Array{Int, 1}=Int[], mesh::Int = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ=  wannier_bands(wannier_file, cell_map_file, [x_mesh/mesh, y_mesh/mesh, 0], nbands)
            for band in 1:nbands
                if band ∉ exclude_bands
                    ϵ_band = ϵ[band]
                    WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]=WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]+histogram_width*(1/mesh)^2
            
                end
        
            end
    
        end
    end

    return WannierDOS

end

function density_of_states_wannier(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Int; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ϵ=  wannier_bands(HWannier, cell_map, [x_mesh/mesh, y_mesh/mesh, 0], nbands)
            for band in 1:nbands
                if band ∉ exclude_bands

                    ϵ_band = ϵ[band]
                    WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]=WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]+histogram_width*(1/mesh)^2
            
                end
        
            end
    
        end
    end

    return WannierDOS

end

function density_of_states_montecarlo(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Int; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)

    xmesh = rand(mesh)
    ymesh = rand(mesh)

    for x_iter in xmesh
        for y_iter in ymesh
            
            ϵ=  wannier_bands(HWannier, cell_map, [x_iter, y_iter, 0], nbands)
            for band in 1:nbands
                if band ∉ exclude_bands

                    ϵ_band = ϵ[band]
                    WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]=WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]+histogram_width*(1/mesh)^2
            
                end
        
            end
    
        end
    end

    return WannierDOS

end

function density_of_states_montecarlo_3d(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Int; exclude_bands = Int[], mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)

    WannierDOS=np.zeros(histogram_width*energy_range)

    xmesh = rand(mesh)
    ymesh = rand(mesh)
    zmesh = rand(mesh)

    for x_iter in xmesh
        for y_iter in ymesh
            for z_iter in zmesh
            
                ϵ=  wannier_bands(HWannier, cell_map, [x_iter, y_iter, 0], nbands)
                for band in 1:nbands
                    if band ∉ exclude_bands

                        ϵ_band = ϵ[band]
                        WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]=WannierDOS[round(Int, histogram_width*(ϵ_band+offset))]+histogram_width*(1/mesh)^2
            
                    end
                
                end
            end
        end
    end

    return WannierDOS

end




function find_chemical_potential(wannier_file::String, cell_map_file::String; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)
    
    doss = density_of_states_wannier(wannier_file, cell_map_file, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    totalstates = []
    for i in 1:length(doss)
        push!(totalstates, [i/histogram_width-offset, sum(doss[1:i]*1/histogram_width)])
    end
    xenergies = []
    yoccupations = []
    for i in 1:length(doss)
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end

    return xenergies, yoccupations

end

function find_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}; mesh::Int=100, histogram_width::Int=100, energy_range::Real=10, offset::Real=0)
    
    doss = density_of_states_wannier(HWannier, cell_map, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    totalstates = []
    for i in 1:length(doss)
        push!(totalstates, [i/histogram_width-offset, sum(doss[1:i]*1/histogram_width)])
    end
    xenergies = []
    yoccupations = []
    for i in 1:length(doss)
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end

    return xenergies, yoccupations

end

function find_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, nbands::Int; mesh::Real = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)
    
    doss = density_of_states_wannier(HWannier, cell_map, nbands, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    totalstates = []
    for i in 1:length(doss)
        push!(totalstates, [i/histogram_width-offset, sum(doss[1:i]*1/histogram_width)])
    end
    xenergies = []
    yoccupations = []
    for i in 1:length(doss)
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end

    return xenergies, yoccupations

end


"""
Returns an array of occupations. Method to find chemical potential at finite temperature. 
"""
function finite_temperature_chemical_potential(HWannier::Array{Float64, 3}, cell_map::Array{Float64, 2}, T::Real; mesh::Real = 100, histogram_width::Real = 100, energy_range::Real = 10, offset::Real = 0)
    
    doss = density_of_states_wannier(HWannier, cell_map, mesh=mesh, histogram_width=histogram_width, energy_range=energy_range, offset=offset )
    occupations_array = Float64[]

    for i in 1:length(doss)
        μ = i/histogram_width-offset
        Fermi = x -> 1/(exp((x-μ)/(kB*T))+1)

        Occupations= Fermi.((1/histogram_width-offset):1/histogram_width:(length(doss)/histogram_width-offset))
        push!(occupations_array, sum(Occupations.*doss)*1/histogram_width)
    end
    
    return collect((1/histogram_width-offset):1/histogram_width:(length(doss)/histogram_width-offset)), occupations_array

end


function find_num_phonons(cell_map::String, phononOmegaSq::String; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 2)
    
    doss = phonon_density_of_states(cell_map, phononOmegaSq; mesh=mesh, histogram_width=histogram_width, energy_range=energy_range)
    totalstates = []
    for i in 1:length(doss)
        push!(totalstates, [i/histogram_width, sum(doss[1:i]*1/histogram_width)])
    end

    xenergies = []
    yoccupations = []
    for i in 1:length(doss)
        push!(xenergies, totalstates[i][1])
        push!(yoccupations, totalstates[i][2])
    end

    return xenergies, yoccupations

end

function density_of_states_wannier(wannier_file_up::String, cell_map_file_up::String, wannier_file_dn::String, cell_map_file_dn::String,; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 10, offset::Real = 0)
   
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

#Next, functions for the calculation of the phonon density of states

"Returns the phonon density of states (phonons per unit energy per unit cell)"
function phonon_density_of_states(cell_map::String, phononOmegaSq::String; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 2)

    PhononDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ωs=  phonon_dispersion(cell_map, phononOmegaSq, [x_mesh/mesh, y_mesh/mesh, 0])
            for ω in ωs
                if ω>0
                    PhononDOS[round(Int, histogram_width*ω)+1]=PhononDOS[round(Int, histogram_width*ω)+1]+histogram_width*(1/mesh)^2
                end
            end
        end
    end

    return PhononDOS


end