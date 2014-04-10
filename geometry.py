import numpy as np
import numpy.linalg as la
import itertools as it

from numpy.random import rand

class UnitCell(object):
  def __init__(self, size, sites): # sites is a list of tuples
    self.size = np.array(size) # unit cell shape
    self.dim = len(self.size.shape)
    self.sites = []
    self.names = []
    for site in sites:
      self.sites.append(np.array(site[0])) # coordination
      self.names.append(sites[1])
    self.nsites = len(self.sites)

class FSuperCell(object):
  def __init__(self, unitcell, size, sites):
    self.unitcell = unitcell
    self.dim = unitcell.dim
    self.csize = np.array(size)
    self.size = self.csize * unitcell.size
    self.ncells = np.product(self.csize)
    self.nsites = unitcell.nsites * self.ncells

    self.sites = []
    self.names = []
    self.unitcell_list = []

    for p in it.product(*tuple([range(a) for a in self.csize])):
      self.unitcell_list.append(np.array(p))
      for i in len(unitcell.sites):
        self.sites.append(unitcell.size * np.array(p) + unitcell.sites[i])
        self.names.append(unitcell.names[i])

    self.fragments = None
    self.Vloc = None
    self.Delta = None

  def set_fragments(self, frag):
    sites = []
    for f in self.frag:
      sites += f.sites
    if sorted(sites) != sorted(self.sites):
      raise Exception("fragment definition incompatible with supercell")
    self.fragments = frag

class FLattice(object):
  def __init__(self, size, sc, bc):
    self.supercell = sc
    self.dim = sc.dim
    if bc == "open":
      self.bc = 0
    elif bc == "pbc":
      self.bc = 1
    elif bc == "apbc":
      self.bc = -1
    else:
      self.bc = bc

    self.scsize = np.array(size)
    self.size = self.scsize * sc.size
    self.nscells = np.product(self.csize)
    self.nsites = sc.nsites * self.nscells

    self.sites = []
    self.names = []
    self.supercell_list = []
    for p in it.product(*tuple([range(a) for a in self.scsize])):
      self.supercell_list.append(np.array(p))
      for i in len(sc.sites):
        self.sites.append(sc.size * np.array(p) + sc.sites[i])
        self.names.append(sc.names[i])
   
    self.h0 = None
    self.h0_kspace = None
  
  def build_h0(self, Ham):
    self.h0 = Ham.build_h0(self)

  def FFTtoK(self, A):
    # currently only for pbc
    assert(self.bc == 1)
    B = A.reshpae(tuple(self.scsize) + A.shape[-2:])
    return np.fft.fftn(B, axes = range(self.dim)).reshape(A.shape)

  def FFTtoT(self, A):
    assert(self.bc == 1)
    B = A.reshape(tuple(self.scsize) + A.shape[-2:])
    C = np.fft.ifftn(B, axes = range(self.dim)).reshape(A.shape)
    if np.allclose(C.imag, 0.):
      return C.real
    else:
      return C

  def get_kpoints(self):
    assert(self.bc == 1)
    self.kpoints = [np.fft.fftfreq(self.scsize[d], 1/(2*np.pi)) for d in range(self.dim)]

  def expand(self, A):
    # expand reduced matrices, eg. Hopping matrix
    assert(self.bc == 1)
    B = np.zeros((self.nsites, self.nsites))
    scnsites = self.supercell.nsites
    nonzero = []
    for i in range(A.shape[0]):
      if not np.allclose(A[i], 0.):
        nonzero.append(i)
    for i in range(self.nscells):
      for j in nonzero:
        # FIXME transform relative position to matrix indice        
        coupled = self.cidx(np.array(self.cpos(i)) + np.array(self.cpos(j)))
        h[i*nfragsite: (i+1)*nfragsite, coupled*nfragsite: (coupled + 1)*nfragsite] = h_reduced[j]
    return h

  def get_neighbor(self, dis = 1., name = None, sites = None):
    # return neighbors
    # for model systems, distance 1. usually means nearest neighbor
    # two lists are returned, first is within the lattice, second is along boundary
    if name is None:
      name = []
      for n in self.supercell.unitcell.names:
        if not n in name:
          name.append(n)
      
    if sites == None:
      sites = [s for s in range(self.nsites) if self.names[s] in name]
    else:
      for s in sites:
        assert(self.names[s] in name)
    
    neighbor1 = [] # in lattice
    neighbor2 = [] # through boundary

    shifts = [np.array(x) for x in it.product([-1, 0, 1], repeat = self.dim) if x != [0] * self.dim]
    for s1 in sites:
      for s2 in range(self.nsites):
        if self.names[s2] in name and s2 > s1:
          if la.norm(self.sites[s2] - self.sites[s1]) - dis < 1e-5:
            neighbor1.append((s1, s2))
          else:
            for shift in shifts:
              if la.norm(self.sites[s2]-self.sites[s1] - shift * self.size) - dis < 1e-5:
                neighbor2.append((s1, s2))
                break

    return neighbor1, neighbor2
