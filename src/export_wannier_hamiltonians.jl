function write_map_write_h(cell_map::String, cell_weights::String, H::String, kmesh::Array{<:Real, 1}, band_file::String, cell_map_file::String)
    py"""   
    def write_map_write_h_py(cell_map, cell_weights, H, kmesh, band_file, cell_map_file):
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

    py"write_map_write_h_py"(cell_map, cell_weights, H, kmesh, band_file, cell_map_file)
end


function write_momentum(cell_map::String, cell_weights::String, H::String, P::String, kmesh::Array{Int, 1}, momentum_file::String)
    py"""

    def write_map_write_p(cell_map, cell_weights, H, P, kmesh, momentum_file):
        import numpy as np
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
        np.savetxt(momentum_file, Pwannier.reshape(len(iReduced), 3*nBands*nBands ))
 
    """
    
    py"write_map_write_p"(cell_map, cell_weights, H, P,  kmesh, momentum_file)
end

function write_eph_matrix_elements(cell_map::String, cell_weights::String, cell_map_ph::String, cell_map_ph_weights::String, HPh::String, nModes::Int, qmesh::Array{Int, 1})
    py"""
    def write_eph(cell_map, cell_weights, cell_map_ph, cell_map_ph_weights, HPh, nModes, qmesh):
        import numpy as np
        cellMap = np.loadtxt(cell_map)[:,0:3].astype(np.int)
        Wwannier = np.fromfile(cell_weights)
        nCells = cellMap.shape[0]
        nBands = int(np.sqrt(Wwannier.shape[0] / nCells))

        cellMapEph = np.loadtxt(cell_map_ph, usecols=[0,1,2]).astype(int)
        nCellsEph = cellMapEph.shape[0]

        prodPhononSup = np.prod(qmesh)
        phononSupStride = np.array([qmesh[1]*qmesh[2], qmesh[2], 1])

        nAtoms = nModes // 3
        cellWeightsEph = np.fromfile(cell_map_ph_weights).reshape((nCellsEph,nBands,nAtoms)).swapaxes(1,2)
        cellWeightsEph = np.repeat(cellWeightsEph.reshape((nCellsEph,nAtoms,1,nBands)), 3, axis=2) #repeat atom weights for 3 directions
        cellWeightsEph = cellWeightsEph.reshape((nCellsEph,nModes,nBands)) #coombine nAtoms x 3 into single dimension: nModes

        iReducedEph = np.dot(np.mod(cellMapEph, qmesh[None,:]), phononSupStride)
        HePhReduced = np.fromfile(HPh).reshape((prodPhononSup,prodPhononSup,nModes,nBands,nBands)).swapaxes(3,4)
        HePhWannier = cellWeightsEph[:,None,:,:,None] * cellWeightsEph[None,:,:,None,:] * HePhReduced[iReducedEph][:,iReducedEph]
        
        return HePhWannier, cellMapEph

    """
    py"write_eph"(cell_map, cell_weights, cell_map_ph, cell_map_ph_weights, HPh, nModes, qmesh)

end



