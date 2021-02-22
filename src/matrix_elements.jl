using PyCall

function phmatrixelements(k1::Array{T, 1}, k2::Array{R, 1}) where {T<:Number, R<:Number}
    py"""
    import numpy as np
    cellMap = np.loadtxt("wannierDefect.mlwfCellMapUp")[:,0:3].astype(int)
    Wwannier = np.fromfile("wannierDefect.mlwfCellWeightsUp")
    nCells = cellMap.shape[0]
    nBands = int(np.sqrt(Wwannier.shape[0] / nCells))
    Wwannier = Wwannier.reshape((nCells,nBands,nBands)).swapaxes(1,2)
    #--- Get k-point folding from totalE.out:
    kfold = np.array([6, 6, 1])
    kfoldProd = np.prod(kfold)
    kStride = np.array([kfold[1]*kfold[2], kfold[2], 1])
    #--- Read reduced Wannier Hamiltonian, momenta and expand them:
    Hreduced = np.fromfile("wannierDefect.mlwfHUp").reshape((kfoldProd,nBands,nBands)).swapaxes(1,2)
    iReduced = np.dot(np.mod(cellMap, kfold[None,:]), kStride)
    Hwannier = Wwannier * Hreduced[iReduced]

    #Read phonon dispersion relation:
    cellMapPh = np.loadtxt('BN33BC.phononCellMap', usecols=[0,1,2]).astype(int)
    nCellsPh = cellMapPh.shape[0]
    omegaSqR = np.fromfile('BN33BC.phononOmegaSq') #just a list of numbers
    nModes = int(np.sqrt(omegaSqR.shape[0] // nCellsPh))
    omegaSqR = omegaSqR.reshape((nCellsPh, nModes, nModes)).swapaxes(1,2)

    #Read e-ph matrix elements
    cellMapEph = np.loadtxt('wannierDefect.mlwfCellMapPhUp', usecols=[0,1,2]).astype(int)
    nCellsEph = cellMapEph.shape[0]
    #--- Get phonon supercell from phonon.out:
    phononSup = np.array([1, 1, 1])
    prodPhononSup = np.prod(phononSup)
    phononSupStride = np.array([phononSup[1]*phononSup[2], phononSup[2], 1])
    #--- Read e-ph cell weights:
    nAtoms = nModes // 3
    cellWeightsEph = np.fromfile("wannierDefect.mlwfCellWeightsPhUp").reshape((nCellsEph,nBands,nAtoms)).swapaxes(1,2)
    cellWeightsEph = np.repeat(cellWeightsEph.reshape((nCellsEph,nAtoms,1,nBands)), 3, axis=2) #repeat atom weights for 3 directions
    cellWeightsEph = cellWeightsEph.reshape((nCellsEph,nModes,nBands)) #coombine nAtoms x 3 into single dimension: nModes
    #--- Read, reshape and expand e-ph matrix elements:
    iReducedEph = np.dot(np.mod(cellMapEph, phononSup[None,:]), phononSupStride)
    HePhReduced = np.fromfile('wannierDefect.mlwfHePhUp').reshape((prodPhononSup,prodPhononSup,nModes,nBands,nBands)).swapaxes(3,4)
    HePhWannier = cellWeightsEph[:,None,:,:,None] * cellWeightsEph[None,:,:,None,:] * HePhReduced[iReducedEph][:,iReducedEph]

    #Constants / calculation parameters:
    eV = 1/27.2114 #in Hartrees

    #Calculate energies, eigenvectors and velocities for given k
    def calcE(k):
        #Fourier transform to k:
        phase = np.exp((2j*np.pi)*np.dot(k,cellMap.T))
        H = np.tensordot(phase, Hwannier, axes=1)
        #Diagonalize and switch to eigen-basis:
        E,U = np.linalg.eigh(H) #Diagonalize
        return E, U

    #Calculate phonon energies and eigenvectors for given q
    def calcPh(q):
        phase = np.exp((2j*np.pi)*np.tensordot(q,cellMapPh.T, axes=1))
        omegaSq, U = np.linalg.eigh(np.tensordot(phase, omegaSqR, axes=1))
        omegaPh = np.sqrt(np.maximum(omegaSq, 0.))
        return omegaPh, U

    #Calculate e-ph matrix elements, along with ph and e energies, and e velocities
    def calcEph(k1, k2):
        #Electrons:
        E1, U1= calcE(k1)
        E2, U2 = calcE(k2)
        #Phonons for all pairs pf k1 - k2:
        omegaPh, Uph = calcPh(k1[:,None,:] - k2[None,:,:])
        #E-ph matrix elements for all pairs of k1 - k2:
        phase1 = np.exp((2j*np.pi)*np.dot(k1,cellMapEph.T))
        phase2 = np.exp((2j*np.pi)*np.dot(k2,cellMapEph.T))
        normFac = np.sqrt(0.5/np.maximum(omegaPh,1e-6))
        g = np.einsum('Kbd,kKycb->kKycd', U2, #Rotate to electron 2 eigenbasis
            np.einsum('kac,kKyab->kKycb', U1.conj(), #Rotate to electron 1 eigenbasis
            np.einsum('kKxy,kKxab->kKyab', Uph, #Rotate to phonon eigenbasis
            np.einsum('KR,kRxab->kKxab', phase2, #Fourier transform from r2 -> k2
            np.einsum('kr,rRxab->kRxab', phase1.conj(), #Fourier transform from r1 -> k1
            HePhWannier))))) * normFac[...,None,None] #Phonon amplitude factor
        return g#, omegaPh
    """
    py"calcEph"(k1, k2) 

end