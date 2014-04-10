import numpy as np
import numpy.linalg as la
import itertools as it

class FHamHubbard(object):
  def __init__(self, U, t = 1.):
    self.U = U
    self.t = t

  def build_h0(self, lattice):
    bc = lattice.bc
    if bc == 1 or bc == -1: # PBC or APBC, then Hamiltonian is reduced
      nsc = lattice.nscells
      ncellsites = lattice.supercell.nsites
      H0 = np.zeros((nsc, ncellsites, ncellsites))
      pairs = lattice.get_neighbor(sites = range(ncellsites))
      for nn in pairs[0]:
        H0[nn[1] / ncellsites, nn[0], nn[1] % ncellsites] = self.t
      for nn in pairs[1]:
        H0[nn[1] / ncellsites, nn[0], nn[1] % ncellsites] = self.t * bc
    elif bc == 0:
      nsites = lattice.nsites
      H0 = np.zeros((nsites, nsites)) # only have real space representation
      pairs = lattice.get_neighbor()
      for nn in pairs[0]:
        H0[nn] = self.t
    else:
      raise Exception("Unsupported boundary condition")
    
    return H0

  def get_Int2e(self):
    # this should be a general function for all types of Hamiltonians
    # it returns a tuple (Int2e, U)
    return (None, self.U)
