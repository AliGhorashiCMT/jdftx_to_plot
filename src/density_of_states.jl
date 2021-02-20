using Plots
using PyCall
using LinearAlgebra 
using Distances

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