# This file is part of the lattice-dmet program. lattice-dmet is free
# software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software
# Foundation, version 3.
# 
# lattice-dmet is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
# 
# Authors:
#    Gerald Knizia, 2012

#from numpy import *
#from helpers import *
import numpy as np
import scipy.linalg as la
import itertools as it

# TODO: consider replacing UnitCell by just Cell

def is_square(M):
   return len(M.shape) == 2 and M.shape[0] == M.shape[1]

def MakeIntegerLattice(Size):
   """return an array of form np.product(Size) x len(Size),
   representing all integers (ijk..) with
      0 <= i < Size[0],
      0 <= j < Size[1],
      ...
   """
   iTs = []
   for kji in it.product(*[range(o) for o in Size][::-1]):
      ijk = kji[::-1] # i goes fastest, then j, then k
      iTs.append(ijk)
   return np.array(iTs)

class FSiteList(object):
   def __init__(self, SitesOrSiteTypes, Xyzs = None):
      """create a list of sites. This corresponds to (a part of) a basis set.
      Args:
         - SiteTypes: List of user defined objects (e.g., strings or own classes),
           identifying different types of sites.
           --OR--
           Sequence of tuples (SiteType,Xyz) with SiteType as above.
         - Xyzs: len(SiteTypes) x nDimR matrix defining the positions of the
           sites in real space. The real-space dimension (xyz dimension) nDimR
           can be chosen freely, but must be equal for all sites. Will usually
           be equal to the lattice dimension, but need non necessarily be.
      """
      if ( Xyzs is None ):
         # given a sequence of tuple -- unpack the input.
         Xyzs = np.array([o[1] for o in SitesOrSiteTypes])
         SiteTypes = [o[0] for o in SitesOrSiteTypes]
      else:
         SiteTypes = SitesOrSiteTypes

      self.SiteTypes = SiteTypes
      self.Types = self.SiteTypes
      if not isinstance(Xyzs, np.ndarray):
         self.Xyzs = np.array(Xyzs) # copy into array format
      else:
         self.Xyzs = Xyzs # reference the original array.
      assert(self.Xyzs.shape[0] == len(self.SiteTypes))

   def __len__(self):
      return len(self.SiteTypes)
   def __getitem__(self, i):
      if isinstance(i, slice):
         return FSiteList(self.SiteTypes[i], self.Xyzs[i])
      else:
         return (self.SiteTypes[i], self.Xyzs[i])
   def __iter__(self):
      for (Type,Xyz) in zip(self.SiteTypes, self.Xyzs):
         yield (Type,Xyz)
   def __add__(self, dXyz):
      """translate the sites by the given vector; return new sites."""
      return FSiteList(self.SiteTypes, self.Xyzs + dXyz)
   def __str__(self):
      L = []
      L.append("   Index " + (len(self.SiteTypes)*"%6s") % tuple(range(len(self.SiteTypes))))
      L.append("   Type  " + (len(self.SiteTypes)*"%6s") % tuple(self.SiteTypes))
      for ixyz in range(self.Xyzs.shape[1]):
         L.append("   Pos/%s " % "xyz"[ixyz] + (len(self.SiteTypes)*"%6s") % tuple(self.Xyzs[:,ixyz]))
      return "\n".join(L)

   def Repeat(self, Size, LatticeVectors):
      """ repeat the sites of *this Size[i] times in direction
      LatticeVectors[i,:]. Returns (iTs, Ts, NewSiteList), where iTs is an
      prod(size) x LatticeVectors.shape[0] integer array and Ts =
      dot(iTs,LatticeVectors)."""

      # make lists of lattice translations for the repetitions
      iTs = MakeIntegerLattice(Size)
      Ts = np.dot(iTs, LatticeVectors.T)

      #print ":Repeat  Size = %s\niTs = \n%s\nVecs =\n%s\nLattiveVectors:\n%s" % (Size,iTs.T,Ts.T, LatticeVectors)

      # make list of (SiteType, xyz) of all sites in the repeated cell.
      # These are simply the unit-cell sites translated by self.Ts
      SiteXyzs = []
      for T in Ts:
         for UcXyz in self.Xyzs:
            SiteXyzs.append(T + UcXyz)
      SiteTypes = len(Ts) * self.SiteTypes
      NewSites = FSiteList(SiteTypes, SiteXyzs)

      ##print "Uc.Repeat %s:\niTs = %s\nTs=%s\n%s" % (Size,iTs,Ts,NewSites)

      return (iTs,Ts,NewSites)


class FLatticeModel(object):
   """abstract base class for a physical lattice model, defining
   both the lattice and the lattice Hamiltonian."""
   def __init__(self, UnitCell, LatticeVectors, ModelParams=None, MaxRangeT=float("inf"),
         SitesAreSpinOrbitals=False, EnergyFactor=1.):
      """Defines the site sof the lattice. You still need to provide
      a function for returning the actual matrix elements.

      Args:
         UnitCell: an list of tuples (SiteType,xyz) defining the sites in the unit
            cell. The type of SiteType is up to you (e.g., a string or an own
            class). xyz must be a numpy array with the same length as the
            lattice dimension.
         LatticeVectors: a D x R matrix, with D the lattice dimension, and R the
            real-space dimension (you could, for example, have two stacked 2d
            lattices in 3d real space, with lattices in xy direction and
            different z components). v[0,:] is the first lattice vector, v[1,:]
            the second, etc. This defines the lattice translations. The unit-
            cell is understood to be repeated in all directions an infinite
            number of times, by adding integer amounts ijk of v[ijk,:] to the
            xyz parameters of the unit cells.

            Note: Lattice translations are indexed with integer vectors (i,j,k),
            to differenciate them from the real-space coordinates (x,y,z).
         ModelParams: If given, defines the names of calculation parameters on
            which this Hamiltonian/Lattice depend (for information and
            consistency purposes)
         MaxRangeT: The range beyond t_ij matrix elements can be considered
            zero (optimization setting for culling at unit-cell level,
            unused atm).
         SitesAreSpinOrbitals: The given unit-cell does contain spin-orbital
            sites instead of spatial orbitals. Note that in this case it is
            requires that all input sites are such that alpha orbitals come
            on even sites and beta orbitals on odd sites!
         EnergyFactor: By default, energies/extensive quantities are normalized
            to the number of unit cells (=1.). This factor can be used, for
            example, to normalize them to the number of sites instead. All
            energies and numbers of electrons are multiplied by this factor.
      """
      assert(isinstance(LatticeVectors, np.ndarray))
      self.nDim = LatticeVectors.shape[1]

      assert(LatticeVectors.shape[1] >= LatticeVectors.shape[0])
      self.LatticeVectors = LatticeVectors
      self.UnitCell = FSiteList(UnitCell)

      assert(self.UnitCell.Xyzs.shape[1] == LatticeVectors.shape[1])

      self.MaxRangeT = MaxRangeT
      self.ModelParams = ModelParams
      self.EnergyFactor = EnergyFactor

      self.SitesAreSpinOrbitals = SitesAreSpinOrbitals
      if self.SitesAreSpinOrbitals:
         # since all sites must come with A and B indices, the total
         # number must be even.
         assert(len(self.UnitCell) % 2 == 0)
         for i in range(0,len(self.UnitCell),2):
            # A and B sites probably should be at the same place.
            assert(self.UnitCell.Xyzs[i] == self.UnitCell.Xyzs[i+1])

      # we'd expect real or complex floats here.
      #self.ScalarType = type(self.GetTij(UnitCell[0],UnitCell[0]))
      #assert(self.ScalarType in (float,complex,np.float64,np.complex128))

   def GetTij(self, SiteI, SiteJ):
      """Get the core Hamiltonian matrix elements t_ij, where
      SiteI = (SiteTypeI, XyzI) and SiteJ = (SiteTypeJ, XyzJ)"""
      raise Exception("FLatticeModel::GetTij must be implemented in derived classes!")

   def GetUi(self, SiteI):
      """return the size of the on-site Hubbard interaction U on site i.
      SiteI = (SiteTypeI, XyzI)"""
      raise Exception("FLatticeModel::GetUi must be implemented in derived classes!")

   def MakeTijMatrix(self, SitesR, SitesC):
      """return a len(SitesR) x len(SitesC) size of the core Hamiltonian
      matrix."""
      #CoreH = np.zeros((len(SitesR), len(SitesC)), self.ScalarType)
      CoreH = np.zeros((len(SitesR), len(SitesC)))
      for (iSiteC, SiteC) in enumerate(SitesC):
         for (iSiteR, SiteR) in enumerate(SitesR):
            CoreH[iSiteR,iSiteC] = self.GetTij(SiteR, SiteC)
      return CoreH

   #def MakeSuperCellSites(self, iSuperCell):
      #"""return the sites of the model, in the format 'list of (SiteType,Xyz)',
      #for a given super-cell id. The super-cell id is an array of integers
      #with the same dimensionality as the latice."""
      #assert(iSuperCell.shape[1] == self.nDim and iSuperCell.dtype == int)
      #nSuperCells = 1
      #if (len(iSuperCell.shape) == 2):
         #nSuperCells = iSuperCell.shape[0]
      #UnitXyz = np.array([o[1] for o in self.UnitCell])
      #UnitSiteTypes = [o[0] for o in self.UnitCell]

      #LatticeTranslations = np.dot(iSuperCell, self.LatticeVectors)

      #SuperCellXyz = []
      #for iSc in range(nSuperCells):
         #lt = LatticeTranslations[iSc,:]
         #for iUnitSite in range(len(UnitSiteTypes)):
            #SuperCellXyz.append(UnitXyz[iUnitSite,:] + lt)
      #return FSiteList(nSuperCells*UnitSiteTypes, SuperCellXyz)


# what follows now are some example classes for simple lattice models which
# demonstrate how the lattice class is intented to be used. For making new
# classes, you probably want to define your own files.

class FHubbardModel1d(FLatticeModel):
   """define a cluster-impurity version of the infinite 1d Hubbard model"""
   def __init__(self, t, U):
      # we give all sites the same type ('X'), and give them a physical
      # location which is simply their lattice index.
      UnitCell = []
      UnitCell.append(("X", np.array([0])))

      self.t = t
      self.U = U

      # The elementary lattice translation is a 1 x 1 matrix: Translating
      # the lattice by one cluster size leaves it invariant.
      LatticeVectors = np.array([[1]])
      FLatticeModel.__init__(self, UnitCell, LatticeVectors, ["U", "t"], MaxRangeT=1)

   def GetTij(self, (SiteTypeI,XyzI), (SiteTypeJ,XyzJ)):
      # return -t/conj(-t) for nearest neighbors in positive/negative x direction.
      dXyz = XyzJ - XyzI
      if dXyz[0] == 1: return -self.t
      if dXyz[0] == -1: return np.conj(-self.t)
      return 0.

   def GetUi(self, (SiteTypeI,XyzI)):
      return self.U

class FHubbardModel2dSquare(FLatticeModel):
   """define the infinite 2d Hubbard model on a square lattice"""
   def __init__(self, ClusterSize, t, U):
      UnitCell = []
      UnitCell.append(("X", np.array([0,0])))

      self.t = t
      self.U = U
      assert(np.isreal(self.t))

      # The elementary lattice translation is a 2 x 2 matrix: Translating
      # the lattice by one cluster size leaves it invariant.
      LatticeVectors = np.array([[1, 0], [0, 1]])
      FLatticeModel.__init__(self, UnitCell, LatticeVectors, ["U", "t"], MaxRangeT=1)

   def GetTij(self, (SiteTypeI,XyzI), (SiteTypeJ,XyzJ)):
      dXyz = XyzJ - XyzI
      if sum(abs(dXyz)) == 1: return self.t
      return 0.

   def GetUi(self, (SiteTypeI,XyzI)):
      return self.U

class FHubbardModel2dTilted(FLatticeModel):
   """The infinite 2d Hubbard model on a 45deg tilted square lattice
   (has less degeneracy than the direct one)"""
   def __init__(self, t, U):
      UnitCell = []
      UnitCell.append( ("U1", np.array([0,0])) )
      UnitCell.append( ("U2", np.array([1,0])) )
      self.t = t
      self.U = U

      LatticeVectors = np.zeros((2,2),int)
      LatticeVectors[0,:] = [1,1]
      LatticeVectors[1,:] = [1,-1]

      # FIXME: MaxRangeT handling is wrong for tilted lattices
      # in super-cell code.
      FLatticeModel.__init__(self, UnitCell, LatticeVectors, ["U", "t"], MaxRangeT=1.5, EnergyFactor=.5)
   def GetTij(self, (SiteTypeI,XyzI), (SiteTypeJ,XyzJ)):
      dXyz = XyzJ - XyzI
      if sum(abs(dXyz)) == 1: return -self.t
      return 0.
   def GetUi(self, (SiteTypeI,XyzI)):
      return self.U



class FTestModel2d(FLatticeModel):
   """a model with non-trivial non-next neighbor interactions; used for testing
   super-cell/symmetry adaption code"""
   def __init__(self, ClusterSize, Params):
      UnitCell = []
      for i,j in it.product(*[range(o) for o in ClusterSize]):
         UnitCell.append(("X", np.array([i,j])))
      self.t = Params['t']
      self.U = Params['U']
      assert(np.isreal(self.t))
      # The elimentary lattice translation is a 1 x 1 matrix: Translating
      # the lattice by one cluster size leaves it invariant.
      LatticeVectors = np.array([[ClusterSize[0], 0], [0, ClusterSize[1]]])
      #FLatticeModel.__init__(self, UnitCell, LatticeVectors, ["U", "t"], MaxRangeT=1)
      FLatticeModel.__init__(self, UnitCell, LatticeVectors, ["U", "t"], MaxRangeT=5)

   def GetTij(self, (SiteTypeI,XyzI), (SiteTypeJ,XyzJ)):
      dXyz = XyzJ - XyzI
      r = np.dot(dXyz,dXyz)
      if ( r > 4**2 ): return 0.
      return np.exp(-r)

   def GetUi(self, (SiteTypeI,XyzI)):
      #return self.U
      return self.U + .2*sum(abs(XyzI))





def _TestLattices():
   L = 12
   nImp = 2
   Hub1d = FHubbardModel1d(nImp, t=1., U=4.)

   Sites = Hub1d.MakeSuperCellSites(np.arange(L/nImp).reshape((L/nImp,1)))
   CoreH = Hub1d.MakeTijMatrix(Sites,Sites)
   print CoreH


if __name__ == "__main__":
   _TestLattices()




