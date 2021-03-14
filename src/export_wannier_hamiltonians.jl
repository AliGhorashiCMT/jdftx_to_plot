const np = PyNULL()
const interpol = PyNULL()
const pyintegrate = PyNULL()
export np, interpol, pyintegrate

function __init__()

    copy!(np, pyimport("numpy"))
    copy!(interpol, pyimport("scipy.interpolate"))
    copy!(pyintegrate, pyimport("scipy.integrate"))

    py"""   
    def write_map_write_h(cell_map, cell_weights, H, kmesh, band_file, cell_map_file):
        import numpy as np
        cellMapUp = np.loadtxt(cell_map)[:,0:3].astype(np.int)
        WwannierUp = np.fromfile(cell_weights)
        nCellsUp = cellMapUp.shape[0]
        nBandsUp = int(np.sqrt(WwannierUp.shape[0] / nCellsUp))
        WwannierUp = WwannierUp.reshape((nCellsUp,nBandsUp,nBandsUp)).swapaxes(1,2)
        kfold=np.array([kmesh[0], kmesh[1], kmesh[2]])
        kfoldProd = np.prod(kfold)
        kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])

        HreducedUp = np.fromfile(H).reshape((kfoldProd,nBandsUp,nBandsUp)).swapaxes(1,2)
        iReducedUp = np.dot(np.mod(cellMapUp, kfold[None,:]), kStride)
        HwannierUp = WwannierUp * HreducedUp[iReducedUp]

        np.savetxt(band_file, HwannierUp.reshape(len(iReducedUp), nBandsUp*nBandsUp ))
        np.savetxt(cell_map_file, cellMapUp)
    """

    py"""

    def write_map_write_p(cell_map, cell_weights, H, P, kmesh, band_file, cell_map_file):

        cellMap = np.loadtxt(cell_map)[:,0:3].astype(np.int)
        Wwannier = np.fromfile(cell_weights)
        nCells = cellMap.shape[0]
        nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
        Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)

        kfold=np.array([kmesh[0], kmesh[1], kmesh[2]])
        kfoldProd = np.prod(kfold)
        kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
        Hreduced = np.fromfile(H).reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)

        Preduced = np.fromfile(P).reshape((kfoldProd,3,nBands,nBands)).swapaxes(2,3)
        iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
        Hwannier = Wwannier * Hreduced[iReduced]
        Pwannier = Wwannier[:,None] * Preduced[iReduced]
        np.savetxt(momentum_file, Pwannier.reshape(len(iReducedUp), 3*nBands*nBands ))
 
    """
end

function write_momentum(cell_map::String, cell_weights::String, H::String, P::String, kmesh::Array{Int, 1}, momentum_file::String)
    py"write_map_write_p"(cell_map, cell_weights, H, P,  kmesh, band_file, map_file)
end

function write_map_write_h(cell_map::String, cell_weights::String, H::String, kmesh::Array{Int, 1}, band_file::String, map_file::String)
    py"write_map_write_h"(cell_map, cell_weights, H, kmesh, band_file, map_file)
end



