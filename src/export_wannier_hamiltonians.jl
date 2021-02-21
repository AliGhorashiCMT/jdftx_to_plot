using PyCall

function __init__()
    py"""   
    def write_map_write_h(cell_map, cell_weights, H, kmesh)
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

        np.savetxt("FullBandsDefect.txt", HwannierUp.reshape(43, 1 ))
        np.savetxt("FullBandsCellMapDefect.txt", cellMapUp)
    """
end
function oneoneone(x)
    print("julia\n")
    return py"oneoneone"(x)
end