function dft_graphene_dos_per_area(;kwargs...)
    a = 1.42*sqrt(3)*1/bohrtoangstrom
    graphene_lattice = lattice([a -a/2 0; 0 a*sqrt(3)/2 0; 0 0 20])
    graphene_ucell_area = unit_cell_area(graphene_lattice)
    DOS_DATA_PATH = joinpath(@__DIR__, "../../data/graphene.in.dos")
    print(DOS_DATA_PATH)
    plot(np.loadtxt(DOS_DATA_PATH)[:, 1]*27.2, np.loadtxt(DOS_DATA_PATH)[:, 2]/27.2/graphene_ucell_area, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
end

function dft_graphene_phonon_dispersion(q::Array{<:Real, 1})

    cell_map_path = joinpath(@__DIR__, "../../data/graphene.in.phononCellMap")
    phonon_omegasq_path = joinpath(@__DIR__, "../../data/graphene.in.phononOmegaSq")
    graphene_force_matrix, graphene_cell_map = phonon_force_matrix(cell_map_path, phonon_omegasq_path)
    return phonon_dispersion(graphene_force_matrix, graphene_cell_map, q)

end

function dft_graphene_wannier_dispersion()
    bands = zeros(8, 300)
    bands_dir = joinpath(@__DIR__, "../../data/wannierbands.txt")
    map_dir = joinpath(@__DIR__, "../../data/wanniercellmap.txt")
    for i in 1:100
        println(i);bands[:, i] = wannier_bands(bands_dir, map_dir, [0, 0.5*i/100, 0], 8)
    end
    for i in 1:100
        println(i+100);bands[:, i+100] = wannier_bands(bands_dir, map_dir, [2/3*i/100, 0.5-(0.5+1/3)*i/100, 0], 8)
    end
    for i in 1:100
        println(i+200);bands[:, i+200] = wannier_bands(bands_dir, map_dir, [2/3-2/3*i/100, -1/3+1/3*i/100, 0], 8)
    end
    return bands
end

function graphene_eph_matrix_elements(k1::Array{<:Real, 1}, k2::Array{<:Real, 1})

    bands_dir = joinpath(@__DIR__, "../../data/wannierbands.txt")
    map_dir = joinpath(@__DIR__, "../../data/wanniercellmap.txt")

    cell_map_dir = joinpath(@__DIR__, "../../data/wannier.graphene.in.mlwfCellMap")
    cell_weights_dir = joinpath(@__DIR__, "../../data/wannier.graphene.in.mlwfCellWeights")

    cell_mapph_dir = joinpath(@__DIR__, "../../data/wannier.graphene.in.mlwfCellMapPh")
    cell_weightsph_dir = joinpath(@__DIR__, "../../data/wannier.graphene.in.mlwfCellWeightsPh")
    HePh_dir = joinpath(@__DIR__, "../../data/wannier.graphene.in.mlwfHePh")

    phonon_cellmap_dir = joinpath(@__DIR__, "../../data/graphene.in.phononCellMap")
    phonon_omegasq_dir = joinpath(@__DIR__, "../../data/graphene.in.PhononOmegaSq")

    HePhWannier, cellMapEph=write_eph_matrix_elements(cell_map_dir, cell_weights_dir, cell_mapph_dir, cell_weightsph_dir, HePh_dir, 6, [2, 2, 1])
    forcemat, mapph = phonon_force_matrix(phonon_cellmap_dir, phonon_omegasq_dir)

    return abs.(eph_matrix_elements(HePhWannier, cellMapEph, forcemat, mapph, bands_dir, map_dir, k1, k2, 8)[4, 5, :])

end