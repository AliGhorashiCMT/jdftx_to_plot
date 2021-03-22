function dft_graphene_dos(;kwargs...)
    a = 1.42*sqrt(3)*1/bohrtoangstrom
    graphene_lattice = lattice([a -a/2 0; 0 a*sqrt(3)/2 0; 0 0 20])
    DOS_DATA_PATH = joinpath(@__DIR__, "../../data/graphene.in.dos")
    print(DOS_DATA_PATH)
    plot(np.loadtxt(DOS_DATA_PATH)[:, 1]*27.2, np.loadtxt(DOS_DATA_PATH)[:, 2]/27.2, linewidth=4, size=(800, 400), xlims = (-2,-0.5), ylims = (0,500/27.2), label="Spin Unpolarized"; kwargs...)
end