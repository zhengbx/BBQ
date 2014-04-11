#
# File: geometry.py
# Author: Bo-Xiao Zheng <boxiao@princeton.edu>
#

import numpy as np
import numpy.linalg as la
import itertools as it

class FUnitCell(object):
    def __init__(self, size, sites): # sites is a list of tuples
        self.size = np.array(size) # unit cell shape
        assert(self.size.shape[0] == self.size.shape[1])
        self.dim = self.size.shape[0]
        self.sites = []
        self.names = []
        for site in sites:
            self.sites.append(np.array(site[0])) # coordination
            self.names.append(site[1])
        self.nsites = len(self.sites)

class FSuperCell(object):
    def __init__(self, unitcell, size):
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
            for i in range(len(unitcell.sites)):
                self.sites.append(np.dot(np.array(p), unitcell.size)  + unitcell.sites[i])
                self.names.append(unitcell.names[i])
  
        self.fragments = None
        self.Vloc = None
        self.Delta = None
  
    def set_fragments(self, frag):
        sites = []
        for f in frag:
            sites += f.sites.tolist()
        if sorted(sites) != range(self.nsites):
            raise Exception("fragment definition incompatible with supercell")
        self.fragments = frag

class FLattice(object):
    def __init__(self, size, sc, bc):
        self.supercell = sc
        self.dim = sc.dim
  
        self.scsize = np.array(size)
        self.size = self.scsize * sc.size
        self.nscells = np.product(self.scsize)
        self.nsites = sc.nsites * self.nscells
  
        if bc == "open":
            assert(self.nscells == 1)
            self.bc = 0
        elif bc == "pbc":
            self.bc = 1
        elif bc == "apbc":
            self.bc = -1
            print "Warning: currently not fully supported"
        else:
            raise Exception("Unsupported Boudary condition")
  
        self.sites = []
        self.names = []
        self.supercell_list = []
        for p in it.product(*tuple([range(a) for a in self.scsize])):
            self.supercell_list.append(np.array(p))
            for i in range(len(sc.sites)):
                self.sites.append(np.dot(np.array(p), sc.size)  + sc.sites[i])
                self.names.append(sc.names[i])
     
        self.h0 = None
        self.h0_kspace = None
        self.neighbor1 = None
        self.neighbor2 = None
    
    def set_Hamiltonian(self, Ham):
        self.Ham = Ham
  
    def get_h0(self, kspace = False):
        if kspace:
            if self.h0_kspace is None:
                self.h0_kspace = self.FFTtoK(self.get_h0())
            return self.h0_kspace
        else:
            if self.h0 is None:
                self.h0 = self.Ham.build_h0(self)
            return self.h0
  
    def FFTtoK(self, A):
        # currently only for pbc
        assert(self.bc == 1)
        B = A.reshape(tuple(self.scsize) + A.shape[-2:])
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
        kpoints = [np.fft.fftfreq(self.scsize[d], 1/(2*np.pi)) for d in range(self.dim)]
        return kpoints
  
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
                pos = (self.supercell_list[i] + self.supercell_list[j]) % self.scsize
                idx = self.supercell_list.index(pos)
                B[i*scnsites:(i+1)*scnsites, idx*scnsites:(idx+1)*scnsites] = A[j]
  
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
        
        neighbor1 = []
        neighbor2 = []
  
        shifts = [np.array(x) for x in it.product([-1, 0, 1], repeat = self.dim) if x != [0] * self.dim]
        # first find supercell neighbors
        for s1 in sites:
            for s2 in range(self.nsites):
                if self.names[s2] in name:
                    if abs(la.norm(self.sites[s2] - self.sites[s1]) - dis) < 1e-5:
                        neighbor1.append((s1, s2))
                    else:
                        for shift in shifts:
                            if abs(la.norm(self.sites[s2]-self.sites[s1] - np.dot(shift, self.size)) - dis) < 1e-5:
                                neighbor2.append((s1, s2))
                                break
  
        return neighbor1, neighbor2

def BuildLatticeFromInput(inp_geom):
    unit = FUnitCell(inp_geom.UnitCell["Shape"],
                     inp_geom.UnitCell["Sites"])
    sc = FSuperCell(unit, np.array(inp_geom.ClusterSize))
    from fragments import FFragment
    frags = [FFragment(f["Sites"], f["ImpSolver"], f["Fitting"])
             for f in inp_geom.Fragments]
    sc.set_fragments(frags)

    assert(np.allclose(np.array(inp_geom.LatticeSize) % np.array(inp_geom.ClusterSize), 0.))
    lattice = FLattice(np.array(inp_geom.LatticeSize)/np.array(inp_geom.ClusterSize),
                       sc, inp_geom.BoundaryCondition)
    return lattice

if __name__ == "__main__":
    # build a 2d square lattice
    sites = [(np.array([0., 0.]), "X")]
    shape = np.array([
      [1., 0.],
      [0., 1.],
    ])
    unit = FUnitCell(shape, sites)
    print "UnitCell"
    print  unit.size
    print "Sites:"
    for i in range(len(unit.sites)):
        print unit.names[i], unit.sites[i], "\t",
        if (i+1)%6 == 0:
            print
    print
    print

    sc = FSuperCell(unit, np.array([2, 2]))
    from fragments import FFragment
    frags = [FFragment([0, 1, 2, 3], "Dmrg", "FullRdm")]
    sc.set_fragments(frags)
    print "SuperCell"
    print sc.size
    print "Sites:"
    for i in range(len(sc.sites)):
        print sc.names[i], sc.sites[i], "\t",
        if (i+1)%6 == 0:
            print
    print
    print
    for i, f in enumerate(sc.fragments):
        print "Fragment", i, f
    print

    lattice = FLattice(np.array([4, 4]), sc, "pbc")
    print "Lattice"
    print lattice.size
    print "Sites:"
    for i in range(len(lattice.sites)):
        print lattice.names[i], lattice.sites[i], "\t",
        if (i+1)%6 == 0:
            print
    print

    print
    print "Lattice Functions"
    print "kpoints"
    print lattice.get_kpoints()
    print "nearest neigbors"
    #print lattice.get_neighbor()[0]
    #print lattice.get_neighbor()[1]
    print "Inside lattice", lattice.get_neighbor(sites = [0,1,2,3])[0]
    print "On boundary   ", lattice.get_neighbor(sites = [0,1,2,3])[1]
    print

    print "Test Hamiltonian Class"

    from ham import FHamHubbard
    hub = FHamHubbard(U = 4., t = 1.)
    lattice.set_Hamiltonian(hub)
    print "Core Hamiltonian (real space)"
    print lattice.get_h0()
    print
    print "Core Hamiltonian (k-space)"
    print lattice.get_h0(kspace = True)
    print
    print "Test FFT"
    print la.norm(lattice.FFTtoT(lattice.FFTtoK(lattice.get_h0())) - lattice.get_h0())
