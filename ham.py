#
# File: geometry.py
# Author: Bo-Xiao Zheng <boxiao@princeton.edu>
#

import numpy as np
import numpy.linalg as la
import itertools as it

class FHamHubbard(object):
    def __init__(self, U, t = 1.):
        self.U = U
        self.t = t

    def build_h0(self, lattice):
        bc = lattice.bc
        if bc in [-1, 0, 1]: # no need to distinguish pbc and obc, because obc is sure to have only one supercell
            nsc = lattice.nscells
            ncellsites = lattice.supercell.nsites
            H0 = np.zeros((nsc, ncellsites, ncellsites))
            pairs = lattice.get_neighbor(sites = range(ncellsites))
            for nn in pairs[0]:
                H0[nn[1] / ncellsites, nn[0], nn[1] % ncellsites] = self.t
            for nn in pairs[1]:
                H0[nn[1] / ncellsites, nn[0], nn[1] % ncellsites] = self.t * bc
        else:
            raise Exception("Unsupported boundary condition")
        
        return H0

    def get_Int2e(self):
        # this should be a general function for all types of Hamiltonians
        # it returns a tuple (Int2e, U)
        return (None, self.U)

def ReadFromDumpFile(FileName):
    with open(FileName, "r") as fdump:
        lines = fdump.readlines()
    
    import re
    tokens = []
    
    # First read header
    for n, line in enumerate(lines):
        
        tokens += re.split(r'[ ,=\n]', line)
        if "&END" in tokens:
            break
    lines = lines[n+1:]
    tokens = [t for t in tokens if t != ""]
    for i in range(len(tokens)):
        if tokens[i] == "NORB":
            nsites = int(tokens[i+1])
        elif tokens[i] == "NELEC":
            nelec = int(tokens[i+1])
        elif tokens[i] == "MS2":
            m = int(tokens[i+1])
    
    nelecA = (nelec+m) / 2
    nelecB = (nelec-m) / 2
    
    # then read integrals
    Int1e = np.zeros((nsites, nsites))
    Int2e = np.zeros((nsites, nsites, nsites, nsites))
    CoreE = 0.
    for line in lines:
        tokens = line.split()
        value = float(tokens[0])
        i = int(tokens[1]) - 1
        j = int(tokens[2]) - 1
        k = int(tokens[3]) - 1
        l = int(tokens[4]) - 1
        if i == j == k == l == -1:
            CoreE += value
        elif k == -1 and l == -1:
            # one electron term
            Int1e[i, j] = Int1e[j, i] = value
        else:
            Int2e[i,j,k,l] = Int2e[i,j,l,k] = Int2e[j,i,k,l] = Int2e[j,i,l,k] = Int2e[k,l,i,j] = Int2e[l,k,i,j] = Int2e[k,l,j,i] = Int2e[l,k,j,i] = value

    return FHamQC(nsites = nsites, Int1e = Int1e, Int2e = Int2e, nelecA = nelecA, nelecB = nelecB, CoreE = CoreE)

class FHamQC(object):
    """
    Quantum chemistry Hamiltonian without periodic boundary condition
    """
    def __init__(self, nsites = None, Int1e = None, Int2e = None, nelecA = None, nelecB = None, CoreE = None, DumpFile = None):
        if nsites is not None:
            self.nsites = nsites
            self.Int1e = Int1e
            self.Int2e = Int2e
            self.nelecA = nelecA
            self.nelecB = nelecB
            self.CoreE = CoreE
        elif DumpFile is not None:
            self = ReadFromDump(DumpFile)
        else:
            raise Exception("Unable to initialize Hamiltonian class")
    
    def build_h0(self, lattice):
        """
        A proper quantum chemistry Hamiltonian is symmetric in spin
        """
        assert(lattice.bc == 0 and lattice.nscells == 1 and lattice.supercell.nsites == self.nsites)
        H0 = np.zeros((1, self.nsites, self.nsites))
        H0[0] = self.Int1e

        return H0

    def get_Int2e(self):
        """
        Int2e is assumed to be a 4-dimensional np.array, with spatial orbitals as indices
        It has to be transformed to spin-orbitals if the embedding system doesn't have spin symmetry
        """
        return (self.Int2e, None)

class FHamQCCrystal(object):
    """
    realistic quantum chemistry Hamiltonian on lattice, e.g. molecular crystal, cuprates
    """
    pass


def Hamiltonian(inp_ham):
    if inp_ham.Type == "Hubbard":
        return FHamHubbard(inp_ham.U, t = 1.)
    else:
        raise KeyError('key %s not exists' % inp_ham.Type)


if __name__ == "__main__":
    ReadFromDumpFile("FCIINP_CN_sto3g_rhf")
