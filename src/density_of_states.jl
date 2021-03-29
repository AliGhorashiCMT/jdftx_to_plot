function bandsoverlayedDOS(dosfile::String, band_file::String, num_bands::Int, num_points::Int, energy_range::Tuple{<:Real, <:Real})
    reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
    exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*1/eV;
    exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*1/eV;

    A = plot(exactenergiesdown, color="black", label="", linewidth=2, ylims = collect(energy_range))
    B = plot!(exactenergiesup, color="purple", label="", linewidth=2, ylims = collect(energy_range))

    C = try
            lowerDOS = argmin(abs.(np.loadtxt(dosfile)[:, 1]*1/eV .- energy_range[1]))
            upperDOS = argmin(abs.(np.loadtxt(dosfile)[:, 1]*1/eV .- energy_range[2]))

            plot(np.loadtxt(dosfile_1)[:, 2]*eV, np.loadtxt(dosfile)[:, 1]*1/eV, linewidth=2, ylims = collect(energy_range), xlims = [0, maximum((np.loadtxt(dosfile)[:, 2]*eV)[lowerDOS:upperDOS]) ])
    catch ##
            lowerDOS = argmin(abs.(np.loadtxt(dosfile, skiprows=1)[:, 1]*1/eV .- energy_range[1]))
            upperDOS = argmin(abs.(np.loadtxt(dosfile, skiprows=1)[:, 1]*1/eV .- energy_range[2]))

            plot(np.loadtxt(dosfile, skiprows=1)[:, 2]*eV, np.loadtxt(dosfile, skiprows=1)[:, 1]*1/eV, linewidth=2, ylims = collect(energy_range), xlims = [0, maximum((np.loadtxt(dosfile, skiprows=1)[:, 2]*eV)[lowerDOS:upperDOS]) ])
    end
    plot(B, C, size = (700, 500), legend = false)
end

function bandsoverlayedDOS2(dosfile1::String, dosfile2::String, band_file::String, num_bands::Int, num_points::Int, energy_range::Tuple{<:Real, <:Real})
    reshaped=reshape(read!(band_file, Array{Float64}(undef, num_bands*num_points*2 )),(num_bands, num_points*2));
    exactenergiesup=permutedims(reshaped, [2, 1])[1:num_points, :]*1/eV;
    exactenergiesdown=permutedims(reshaped, [2, 1])[num_points+1:2*num_points, :]*1/eV;

    A = plot(exactenergiesdown, color="black", label="", linewidth=2, ylims = collect(energy_range))
    B = plot!(exactenergiesup, color="purple", label="", linewidth=2, ylims = collect(energy_range))

    ##Load DOS 1
    
    dosdata1 = try 
        np.loadtxt(dosfile1)
    catch 
        np.loadtxt(dosfile1, skiprows=1)
    end
    ##Load DOS 2
    dosdata2 = try 
         np.loadtxt(dosfile2)
    catch 
        np.loadtxt(dosfile2, skiprows=1)
    end

    lowerDOS1 = argmin(abs.(dosdata1[:, 1]*1/eV .- energy_range[1]))
    upperDOS1 = argmin(abs.(dosdata1[:, 1]*1/eV .- energy_range[2]))

    lowerDOS2 = argmin(abs.(dosdata2[:, 1]*1/eV .- energy_range[1]))
    upperDOS2 = argmin(abs.(dosdata2[:, 1]*1/eV .- energy_range[2]))

    max1 = maximum((dosdata1[:, 2]*eV)[lowerDOS1:upperDOS1])
    max2 = maximum((dosdata2[:, 2]*eV)[lowerDOS2:upperDOS2])

    max = maximum([max1, max2])

    C = plot(dosdata1[:, 2]*eV, dosdata1[:, 1]*1/eV, linewidth=2, ylims = collect(energy_range), xlims = [0, max])
    plot!(dosdata2[:, 2]*eV, dosdata2[:, 1]*1/eV, linewidth=2, ylims = collect(energy_range), xlims = [0, max], legend = false)
    
    plot(B, C, size = (700, 500))
end


function density_of_states(dosfile_1::String, dosfile_2::String; kwargs... )
    
    try
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    catch ##
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, np.loadtxt(dosfile_1, skiprows=1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    end

    try
        plot!(np.loadtxt(dosfile_2)[:, 1]*1/eV, np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    catch
        plot!(np.loadtxt(dosfile_2, skiprows=1)[:, 1]*1/eV, np.loadtxt(dosfile_2, skiprows=1)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    end

end

function density_of_states(dosfile_1::String; kwargs...)
    try
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    catch 
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, np.loadtxt(dosfile_1, skiprows=1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    end
end

"""
The typical density of states outputed by JDFTX is per unit cell. However, sometimes it is more relevant to know the 
density of states per unit volume. This is simply equivalent to dividing the conventional DOS by the unit cell 
volume 
"""
function density_of_states_per_area(dosfile_1::String, lattice_vecs::Array{<:Array{<:Real, 1}, 1}; kwargs...)

    ucell_area = unit_cell_area(lattice_vecs)
    try 
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    end

end

function density_of_states_per_area(dosfile_1::String, lattice_vecs::lattice; kwargs...)

    ucell_area = unit_cell_area(lattice_vecs)
    try
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
    end

end


function density_of_states_per_area(dosfile_1::String, dosfile_2::String, lattice_vecs::Array{<:Array{<:Real, 1}, 1}; kwargs... )
    ucell_area = unit_cell_area(lattice_vecs)
    try 
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    end
    
    try
        plot!(np.loadtxt(dosfile_2)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    catch
        plot!(np.loadtxt(dosfile_2, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    end
end

function density_of_states_per_area(dosfile_1::String, dosfile_2::String, lattice_vecs::lattice; kwargs... )
    ucell_area = unit_cell_area(lattice_vecs)
    try 
        plot(np.loadtxt(dosfile_1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    catch
        plot(np.loadtxt(dosfile_1, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_1)[:, 2]*eV, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Up"; kwargs...)
    end

    try
        plot!(np.loadtxt(dosfile_2)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    catch
        plot!(np.loadtxt(dosfile_2, skiprows=1)[:, 1]*1/eV, 1/ucell_area*np.loadtxt(dosfile_2)[:, 2]*eV, linewidth=4,  size=(800, 400), label="Spin Down"; kwargs...)
    end
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

function find_num_phonons(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 2)
    
    doss = phonon_density_of_states(force_matrix, phonon_cell_map; mesh=mesh, histogram_width=histogram_width, energy_range=energy_range)
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


function phonon_density_of_states(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 2)

    PhononDOS=np.zeros(histogram_width*energy_range)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ωs=  phonon_dispersion(force_matrix, phonon_cell_map, [x_mesh/mesh, y_mesh/mesh, 0])
            for ω in ωs
                if ω>0
                    PhononDOS[round(Int, histogram_width*ω)+1]=PhononDOS[round(Int, histogram_width*ω)+1]+histogram_width*(1/mesh)^2
                end
            end
        end
    end

    return PhononDOS
end

function phonon_density_of_states_per_area(force_matrix::Array{<:Real, 3}, phonon_cell_map::Array{<:Real, 2}, lattice_vecs::Array{<:Array{<:Real, 1}, 1}; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 2)

    PhononDOS=np.zeros(histogram_width*energy_range)

    ucell_area = unit_cell_area(lattice_vecs)

    for x_mesh in 1:mesh
        for y_mesh in 1:mesh
            
            ωs=  phonon_dispersion(force_matrix, phonon_cell_map, [x_mesh/mesh, y_mesh/mesh, 0])
            for ω in ωs
                if ω>0
                    PhononDOS[round(Int, histogram_width*ω)+1]=PhononDOS[round(Int, histogram_width*ω)+1]+histogram_width*(1/mesh)^2
                end
            end
        end
    end

    return PhononDOS/ucell_area
end


function phonon_density_of_states_per_area(cell_map::String, phononOmegaSq::String, lattice_vecs::Array{<:Array{<:Real, 1}, 1}; mesh::Int = 100, histogram_width::Int = 100, energy_range::Real = 2)


    ucell_area = unit_cell_area(lattice_vecs)
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

    return PhononDOS/ucell_area
end