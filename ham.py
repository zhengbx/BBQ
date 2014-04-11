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


class FHamQC(object):
    """
    Quantum chemistry Hamiltonian
    """
    def __init__(self, nsites = None, Int1e = None, Int2e = None, nelecA = None, nelecB = None, DumpFile = None):
        if self.nsites is not None:
            self.nsites = nsites
            self.Int1e = Int1e
            self.Int2e = Int2e
            self.nelecA = nelecA
            self.nelecB = nelecB
        elif DumpFile is not None:
            self.nsites, self.Int1e, self.Int2e, self.nelecA, self.nelecB = ReadFromDump(DumpFile)
        else:
            raise Exception("Unable to initialize Hamiltonian class")
    
    def build_h0(self, lattice):
        assert(lattice.bc == 0 and lattice.nscells == 1 and lattice.supercell.nsites == self.nsites)
        H0 = np.zeros((1, self.nsites, self.nsites))
        H0[0] = self.Int1e

        return H0

    def get_Int2e(self):
        return (self.Int2e, None)

def Hamiltonian(inp_ham):
    if inp_ham.Type == "Hubbard":
        return FHamHubbard(inp_ham.U, t = 1.)
    else:
        raise KeyError('key %s not exists' % inp_ham.Type)
