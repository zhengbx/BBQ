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

import numpy as np
from scipy import linalg as la
import itertools as it
from lattice_model import *
from helpers import mdot, is_hermitian, InvertPermutation

def to_ref_type(Op,dtype):
   """cast a complex number to real, if the complex is, for some reason,
   known to actually be real. ndarray.astype does not quite do that, because
   it emits a warning message even for explicit casts (interestingly, that
   is intended behavior)"""
   if (Op.dtype in (complex,np.complex128)) and\
      (dtype in (float,np.float64)):
      assert(np.allclose(Op.imag,0.0))
      return Op.real
   return Op

class FSuperCellParams(object):
   def __init__(self, TotalSize, PhaseShift=None, ClusterSize=None):
      self.TotalSize = TotalSize
      self.PhaseShift = PhaseShift
      self.ClusterSize = ClusterSize

class FSuperCell(object):
   """A super-cell represents a finite cylic model of a periodic lattice
   of unit-cells.
   Notes:
   - self.Size[i] repetitions of the unit-cell in Lattice.T[i] direction
     are explicitly treated; after that, the previous sites are repeated.
   - The super-cell supplies methods for construction operators on the cyclic
     model, and for transforming them between real-space translations (t-space)
     and k-space.
   - The actual calculation basis is defined as follows. Let T be the
     super-cell translation vectors. Then the actual calculation basis is:
         |i>_sc := \sum_{ijk} (PhaseShift[0] * Size[0] * T[:,0])^i *
                              (PhaseShift[1] * Size[1] * T[:,1])^j *
                              (PhaseShift[2] * Size[2] * T[:,2])^k * |i>_uc
     where ijk run over integers, and T  are the unit-cell translations. |i>_uc
     is a physical site in the unit-cell; the cyclic model is, however, defined
     in terms of fixed linear combinations of such sites with the super-cell
     periodic symmetry.
   """
   def __init__(self, Model, TotalSize, PhaseShift=None, ClusterSize=None, OrbType=None):
      """Construct a finite cyclic model, treating N repetitions of
      Model.UnitCell explicitly.
      Args:
         - Model: Lattice model for which this is a super-cell
         - Size: Size[i] is the repretition count in Lattice.T[:,i] direction.
                 Must be a sequence of the same size as the lattice dimension.
         - ClusterSize: If ClusterSize[i] != 0, this many repetitions of
                 Model's unit cell are handled explicitly in real space
                 instead of via translational symmetry. Useful for setting up
                 impurity cluster calculations of lattice models.
                 ClusterSize[i] must be a divisor of Size[i].
         - PhaseShift: Defines if boundary conditions are periodic, anti-periodic,
                 or twisted. This is the phase factor applied when building the
                 fixed linear combinations of physical sites which form the
                 super-cell basis.
         - OrbType: If "UHF", the model will be treated as a spin-orbital model,
                 even if it does not call for it itself (required for UHF mean-
                 fields). That is, the unit cell sites will be doubled, with
                 even sites standing for alpha orbitals, and odd sites for beta
                 orbitals. If None: Use spin orbitals for spin-orbital models,
                 spatial orbitals for spatial orbital models.
      """
      self.Model = Model

      if Model.LatticeVectors.shape[0] != len(TotalSize):
         raise Exception("Lattice dimension and super-cell size don't match!")
      if ClusterSize is not None and len(ClusterSize) != len(TotalSize):
         raise Exception("ClusterSize dimension and super-cell dimension don't match!")

      # size of the super-cell in terms of the physical model's unit-cell
      self.TotalSize = np.array(TotalSize)

      # check if we should expand the lattice model's physical unit cell into a
      # cluster-unit-cell (in the remaining code, the cluster-unit-cell would be
      # the actual unit-cell, and SuperCell.Size is given in terms of the
      # cluster unit-cell)
      if ClusterSize is None:
         # nope.
         self.ClusterSize = np.ones(len(TotalSize),int)
      else:
         # yes.
         self.ClusterSize = np.array(ClusterSize,int)

      # size of super-cell handled via translational symmetry
      # (i.e., size of te super-cell in terms of the cluster unit-cell)
      assert(len(self.ClusterSize) == len(self.TotalSize) and
             np.all(self.ClusterSize >= 0))
      self.Size = self.TotalSize/self.ClusterSize
      self.nDim = len(self.Size)

      assert(np.all(self.Size * self.ClusterSize == self.TotalSize))

      # make lists of translations and sites in the cluster unit-cell...
      d1,d2,self.UnitCell = Model.UnitCell.Repeat(self.ClusterSize, Model.LatticeVectors)
      self.T = self.ClusterSize * Model.LatticeVectors
      # ^- cluster-cell (that's the unit-cell of the super-cell) translation vectors
      #print "SuperCell-T:\n%s" % self.T
      #raise SystemExit

      # if necessary, double the sites to get from spatial orbitals to
      # spin orbitals.
      assert(OrbType in ["UHF", "RHF", None])
      self.SpinOrbs = (OrbType == "UHF") or self.Model.SitesAreSpinOrbitals
      if self.SpinOrbs and not self.Model.SitesAreSpinOrbitals:
         # create interleaved alpha and beta orbitals:
         #    Site0[A], Site0[B], Site1[A], Site1[B], ...
         # -> alpha sites on even indices, beta sites on odd indices.
         SitesA = list(self.UnitCell)
         SitesB = list(self.UnitCell)
         self.UnitCell = FSiteList(list(it.chain(*zip(SitesA,SitesB))))
      else:
         # if the model calls for spin orbitals, it is not okay to explicitly
         # disable them.
         if self.Model.SitesAreSpinOrbitals:
            assert(self.SpinOrbs)

      # ...and the remaining super-cell (i.e., repeating cluster-cells)
      self.iTs, self.Ts, self.Sites = self.UnitCell.Repeat(self.Size, self.T)
      self.iKs = self.iTs

      # maps translation vectors to their indices. should likely
      # not be here..
      self.Fixme_iTsLookup = dict((tuple(t),it) for (it,t) in enumerate(self.iTs))

      # check & store phase factors. Phase factors must be roots of unity.
      if ( PhaseShift is None ):
         self.PhaseShift = np.ones(self.T.shape[1])
         self.dtype = float
      else:
         self.dtype = type(PhaseShift[0])
         if ( self.dtype is int ):
            self.dtype = float
         self.PhaseShift = np.array(PhaseShift, self.dtype)
      assert(np.allclose(np.conj(self.PhaseShift)*self.PhaseShift - 1.0, 0.0))

      # public attributes.
      self.nUnitCells = self.Ts.shape[0]
      self.nClusterSites = np.product(self.ClusterSize)
      self.nUnitCellsTotal = self.nUnitCells * self.nClusterSites
      self.nSitesU = len(self.UnitCell)
      self.nSites = self.nUnitCells * len(self.UnitCell)

      if 0:
         print "Super-cell sites:"
         print self.Sites
         print "!x total number of sites in super-cell: %i" % len(self.Sites)

      self._CalcSuperCellSymmetrizationRange()
      self._PrepareTranslationTrafos()
      self.EnergyFactor = (1.*self.Model.EnergyFactor)/self.nClusterSites

   def _PrepareTranslationTrafos(self):
      """Calculate intermediate data for performing periodic symmetrization
      opererations (t-space vs k-space transformations and matrix reconstruction
      from symmetry-unique elements)."""
      def MakeElementaryTransOp(N, PhaseFactor):
         # set up a right-shift unit-cell translation operator matrix: That is
         # the matrix mapping a sequence like
         #      [x0 x1 x2 x3 x4 x5 x6 x7]
         # into
         #      [x1 x2 x3 x4 x5 x6 x7 phase-factor*x0]
         # where x* incorporates all degrees of freedom not associated
         # with the shift by one unit-cell in the current super-cell direction
         # (that is: the operator acts as
         #        id(sites-in-the-unit-cell) x id(other-directions))
         # for the other degrees of freedom.)
         TransOp = np.eye(N,dtype=type(PhaseFactor))
         TransOp = np.roll(TransOp, -1, axis=0)
         TransOp[:,:1] *= PhaseFactor
         return TransOp

      # make the translation eigenvectors. These effectively define
      # the symmetry adapted basis (i.e., k-space).
      def MakeEvsForDirection(N, PhaseFactor):
         """Args: N: Size of the super-cell in current direction.
         PhaseFactor: PhaseShift for the current direction.

         Returns 2*pi*k, ew[k] (==exp(1j*2*pi*k)), ev[k]"""
         TransOp = MakeElementaryTransOp(N,PhaseFactor)
         TransOp = np.conj(TransOp) # <- FIXME: should that be here?
         ew,ev = la.eig(TransOp)
         assert(np.allclose(np.abs((np.conj(ew)*ew) - 1.), 0.))
         # ^- eigenvalues should all be roots of unity (TransOp is unitary)

         # convert to phase angles and sort by angle. We sort by
         # absolute value first, and then by sign. The absolute value
         # should be associated with the size of the impulse corresponding
         # to the translation, and the sign with the direction.
         phi = np.angle(ew)
         I = np.lexsort((np.sign(phi),np.abs(phi)))
         ew = ew[I]
         ev = ev[:,I]
         phi = phi[I]
         if 0:
            print "sorted phi:", phi
            print "exp(1j*phi) - ew:", la.norm(np.exp(1j*phi) - ew)
            print "TranOp*Ev - Ev*diag(ew): %s" % la.norm(np.dot(TransOp,ev) - np.dot(ev,np.diag(ew)))
         return phi,ew,ev
      if 0:
         MakeEvsForDirection(6,1.)
         MakeEvsForDirection(6,-1.)
         MakeEvsForDirection(6,np.exp(0.12j))
      self.kTrafos = []
      self.kAngles = []
      self.kEws = []
      self.kEvs = []
      for i in range(len(self.Size)):
         phi,ew,ev = MakeEvsForDirection(self.Size[i], self.PhaseShift[i])
         self.kTrafos += [(phi,ew,ev)]
         self.kAngles += [phi]
         self.kEws += [ew]
         self.kEvs += [ev]

      # make a set of permutations and phase factors, for re-constructing
      # the missing rows/columns from a symmetry-unique subset of a full
      # super-cell operator (e.g., for making all rows of CoreH when
      # only the UnitCell.nSites x SuperCell.nSites subset is given)
      #
      # This effectively means calculating the result (permutation and phase)
      # of applying the elementary translation operator i times,
      # for i = 0 ... N-1.
      def MakePermOpsForDirection(N, PhaseFactor):
         TransOp = MakeElementaryTransOp(N,PhaseFactor)
         TransOp = np.conj(TransOp.T) # <- FIXME: should that be here?
         # make an initial permutation: identity. We will be propagating
         # this permutation and storing the result. We actually absorb
         # both the phases (modulus(P))  and the indices (abs(P)) in the
         # same array during the calculation.
         P = (np.arange(N) + 1).astype(TransOp.dtype)
         I = np.zeros((N,N),int)
         PF = np.zeros((N,N),TransOp.dtype)
         for i in range(N):
            I[:,i] = np.rint(np.abs(P) - 1.)
            if 0:
               # FIXME: is this right? this inverts the permutation;
               # I think applies the operator into the opposite direction.
               # Need to check with non-trivial phases!!
               I_ = list(I[:,i])
               I[:,i] = np.array([I_.index(o) for o in range(len(I_))])
            PF[:,i] = P/np.abs(P)
            P = np.dot(TransOp, P)
         return PF,I
      self.tPermsAndPhases = []
      self.tPermsAndPhasesInvI = []
      for i in range(len(self.Size)):
         PF,I = MakePermOpsForDirection(self.Size[i], self.PhaseShift[i])
         self.tPermsAndPhases += [(PF,I)]
         #print "I.shape = ", I.shape
         #print "I[:,2] =     ", I[:,2]
         #print "I[2,:] =     ", I[2,:]
         #print "Inv(I[:,2]) =", InvertPermutation(I[:,2])
         #raise Exception("absicht.")
         InvI = np.array([InvertPermutation(I[:,i]) for i in range(self.Size[i])])
         self.tPermsAndPhasesInvI += [(PF,InvI)]
      pass

   def _MakeUnitCellPhasesForK(self, ijk):
      """return a complex array of size product(self.Size), giving the
      contributions of the represented unit-cells to k-space vectors ijk
      (ordered as self.Ts). ijk is an integer vector, indexing into
      self.kAngles, and has the same bounds as self.Size"""
      X = np.ones(())
      for i in range(self.nDim)[::-1]:
         ev = self.kEvs[i]
         #print "k-ev, dim i = %i: %s" % (i,ev[:,ijk[i]])
         X = np.outer(X, ev[:,ijk[i]])
      #print "X before flatten (should have correct shape!):\n%s" % X
      X = X.reshape((X.size))
      assert(len(X) == self.nUnitCells)
      return X

   def _MakeRealspaceOrb(self, UnitCellOrb, ijk):
      """make a real-space representation of the bloch wave described by the
      UnitCell.nSites x nOrb matrix UnitCellOrb and the k-space
      vector indexed by ijk. ijk is a integer vector, indexing into
      self.kAngles, and has the same bounds as self.Size.

      Output is indexed as T x UnitCell.nSites x nOrb"""
      X = self._MakeUnitCellPhasesForK(ijk)
      NewShape = [self.nSites] + list(UnitCellOrb.shape[1:])
      return np.outer(X, UnitCellOrb).reshape(*NewShape)

   def _MakeUnitCellPhasesForT(self, ijk):
      """returns (PF, I), where I is an integer array indexing into
      self.nUnitCells, and PF is a real or complex array of the same length
      describing phase fectors. Those arrays have the property that the  ijk't
      column of a real-space periodically symmetric matrix M[nUnitsCells x
      nUnitCells] can be obtained as
         M[:,ijk] == PF * M[I,0].
      ijk should be a tuple from self.iTs."""
      I = np.zeros((),int)
      PF = np.ones((),self.PhaseShift.dtype)
      for i in range(self.nDim)[::-1]:
         PFi,Ii = self.tPermsAndPhases[i]
         I = np.add.outer(I*Ii.shape[0], Ii[:,ijk[i]])
         PF = np.outer(PF, PFi[:,ijk[i]])
      return PF.reshape((self.nUnitCells)),\
              I.reshape((self.nUnitCells))
      pass

   def _ConstructTranslatedRows(self, OpSymUq, ijk):
      """given the SuperCell.nSites x UnitCell.nSites matrix OpSymUnique,
      defining a periodically symmetric operator on via its matrix elements on
      the first unit cell with respect to all other unit cells, construct
      the ijk'th row group where ijk is a vector within self.iTs."""
      PF,I = self._MakeUnitCellPhasesForT(ijk)
      U = len(self.UnitCell)
      T = self.nUnitCells
      #print "PF: %s" % PF
      #print "I:  %s" % I
      #assert(OpSymUq.shape == (U,U*T))
      #return (OpSymUq.reshape(((U*U),T))[:,I] * PF).reshape((U,U*T))
      assert(OpSymUq.shape == (T*U,U))
      return (OpSymUq.reshape((T,(U*U)))[I,:].swapaxes(0,1) * PF).swapaxes(0,1).reshape((T*U,U))

   def _MakeUnitCellPhasesForT_Restricted(self, ijkRow, iOpTs):
      # <r + T|op| s> = <r + T + S|op|s + S> for any unit-cell translation S,T
      # we know the op is only non-zero for a limited range of T, which we
      # have given via iOpTs.
      #
      # we have given the ijkRow (lhs T), and are now looking for the ijkCol
      # (right T) for which <r+Tr| op |r+Tc> is non-zero. That means that
      # Tr-Tc must lie within the Tc supplied in iOpTs.

      # FIXME: do this properly.
      if 0:
         I = []
         PF = []
         for ijkOp in (self.iTs[o] for o in iOpTs):
            for iTCol,ijkCol in enumerate(self.iTs):
               if np.all((ijkCol - ijkRow) % self.Size == ijkOp):
                  I.append(iTCol)
                  PF.append(np.product(self.PhaseShift**((ijkCol - ijkRow) / self.Size)))
         return np.array(PF), np.array(I)
      else:
         # well... that's not properly, but atm it's not the main problem in this form.
         I = []
         PF = []
         for ijkOp in (self.iTs[o] for o in iOpTs):
            ijkCol = ijkOp + ijkRow
            PF.append(np.product(self.PhaseShift**((ijkCol) / self.Size)))
            I.append(self.Fixme_iTsLookup[tuple(ijkCol % self.Size)])
         return np.array(PF), np.array(I)

   def _ExpandFull(self, OpSymUq_):
      #print OpSymUq_.shape
      nSitesU, nSitesS = self.nSitesU, self.nSites
      assert(OpSymUq_.shape == (nSitesS,nSitesU) or OpSymUq_.shape == (self.nUnitCells,self.nSitesU,self.nSitesU))
      OpSymUq = OpSymUq_.reshape(nSitesS,nSitesU)
      #print OpSymUq.shape
      OpFull = np.zeros((nSitesS,nSitesS),complex)
      for (i,ijk) in enumerate(self.iTs):
         OpFull[nSitesU*i:nSitesU*(i+1),:] = self._ConstructTranslatedRows(OpSymUq, ijk).T
      return OpFull

   #def ConvertOpToKSpace(self, OpT):
      #"""convert a periodic 1-electron matrix in real space, explicitly defined on
      #sites I x L, to k-space. Return matrix G x I x I"""
      #assert(OpT.shape == (I,L))
      #OpK = (1.+0.j)*zeros((G,I,I))

      #for (ik,k) in enumerate(kGrid.Points):
         #OpK[ik,:,:] = OpT[:,:I]
         #for i in range(1,L/I):
            #phase = exp(1.j*k*i)            # <- phase between cell #0 and cell #i.
            #T = array(range(i*I, (i+1)*I))  # <- sites of the translated cell
            #OpK[ik,:,:] += phase * OpT[:,T]
            #if not IsHermitian(OpK[ik,:,:]):
               #raise Exception("Input real-space Hamiltonian does not correspond to a periodic system!")
      #return OpK

   def _CalcSuperCellSymmetrizationRange(self):
      # minimum array of super-cell iTs over which we need to
      # sum to cover the whole range of Model.MaxRangeT
      # (the latter is in xyz space)
      self.MinSuperCell_iTs = []
      DistXyz =  lambda dxyz: np.sum(dxyz**2,0)**.5
      #DistXyz = lambda dxyz: np.max(np.abs(dxyz),0)
      #DistXyz = lambda dxyz: np.sum(np.abs(dxyz),0)
      # ^- sum of absolute values of coordinate differences in xyz directions
      #    e.g., a diagonal connection would have distance 2 (not 1).
      if self.Model.MaxRangeT < 2:
         #assert(self.Model.MaxRangeT == 1)
         # one super-cell displacement in each direction in the case
         # of nearest neighbor hopping.
         for ijk_ in it.product(*[[-1,0,1] for o in self.Size]):
            ijk = np.array(list(ijk_))
            if DistXyz(ijk) <= self.Model.MaxRangeT + 1e-8:
               self.MinSuperCell_iTs.append(ijk)
      else:
         # this is the more correct thing to do, but actually quite slow for
         # large super-cells, due to the N^2 nearest/farest neighbor search.
         # Should make better algorithms some time, but in Python that is a
         # problem...
         def Dist(T,Crit=np.min):
            # minimum or maximum distance between the sites in the super-cell
            # and the same sites translated by T
            if 0:
               # requires too much memory in big super-cells
               D,N = self.SiteXyzs.shape
               SitesXyzsT = (self.SiteXyzs.T + T).T
               dXyz = (self.SiteXyzs.reshape(D,N,1) - SitesXyzsT.reshape(D,1,N)).reshape(D,N*N)
               return Crit(DistXyz(dXyz))
            else:
               dmin = None
               for v in self.Sites.Xyzs:
                  dXyz = Crit(DistXyz(  ((v+T) - self.Sites.Xyzs) ))
                  if dmin is None: dmin = dXyz
                  else: dmin = Crit([dmin,dXyz])
               return dmin
         def MinDist(T): return Dist(T,np.min)
         def MaxDist(T): return Dist(T,np.max)

         if 0:
            print "MinDist(2,0): %s" % MinDist(np.array([2,0])*self.Size)
            print "MinDist(1,0): %s" % MinDist(np.array([1,0])*self.Size)
            print "MinDist(0,1): %s" % MinDist(np.array([0,1])*self.Size)

         MaxDist1 = MaxDist(0.*self.Size)
         BaseRange = 4*self.Size+np.ones(len(self.Size),int)
         for ijk_ in it.product(*[range(-o/2,o/2+1) for o in BaseRange]):
            # dXyz: change in xyz due to translation by *whole* super-cell ijk.
            ijk = ijk_ * self.Size
            dXyz = np.dot(self.T, ijk)
            if (DistXyz(dXyz) <= self.Model.MaxRangeT + 1e-8 + MaxDist1 and
               True ):
               #MinDist(dXyz) <= Model.MaxRangeT + 1e-8):
               #print "Test: %s  -> dXyz: %s -> MinDist = %s" % (ijk, dXyz, MinDist(dXyz))
               self.MinSuperCell_iTs.append(ijk_)
         #print "self.MinSuperCell_iTs:\n%s" % self.MinSuperCell_iTs
         print "!x done with supercell-displacements"
      #print "self.MinSuperCell_iTs:\n%s" % self.MinSuperCell_iTs

   def MakeTijMatrix(self, SitesR=None, SitesC=None, _SpinOrb=None):
      # differs from Model.MakeTijMatrix by applying the super-cell
      # periodicity and inserting spin-orbitals if the input model
      # is not already a spin-orbital model.
      #if SitesR is None: SitesR = self.UnitCell
      if SitesR is None: SitesR = self.Sites
      if SitesC is None: SitesC = self.Sites

      if _SpinOrb is None and self.SpinOrbs and not self.Model.SitesAreSpinOrbitals:
         # we're supposed to make spin-orbitals, but the input model
         # is not actually a spin-orbital model. So we just ask it to
         # make a Tij matrix for alpha sites, and then replicate that to
         # the beta sites.
         TijA = self.MakeTijMatrix(SitesR[::2], SitesC[::2], _SpinOrb=False)
         Tij = np.zeros((2*TijA.shape[0], 2*TijA.shape[1]), TijA.dtype)
         Tij[ ::2, ::2] = TijA
         Tij[1::2,1::2] = TijA
         return Tij
      else:
         Tij = None
         for ijk in self.MinSuperCell_iTs:
            Phase = np.product(self.PhaseShift**ijk)
            #print "translation: %s  %s" % (np.dot(self.T, self.Size*ijk), ijk)
            SitesR_Trans = SitesR + np.dot(self.T, self.Size*ijk)
            Tij_iT = Phase * self.Model.MakeTijMatrix(SitesR_Trans, SitesC)
            if ( Tij is not None ):
               Tij += Tij_iT
            else:
               Tij = Tij_iT
         assert(Tij is not None)
         return Tij
      pass

   def MakeUiMatrix(self, Sites=None):
      """return the magnitude of the local Hubbard-type two-electron
      Coulomb interactions on the super-cell sites."""
      if Sites is None: Sites = self.UnitCell
      Int2e_U = np.array([self.Model.GetUi(o) for o in Sites])
      assert(Int2e_U.shape == (len(Sites),))
      return Int2e_U

   def _OpShapeForAxis(self, iAxis):
      """return the shape combination (Lower,Size[iAxis],Upper) for representing
      the T or K components of an operator in the direction iAxis."""
      assert(iAxis < len(self.Size))
      Lower = 1
      Upper = self.nSitesU * self.nSitesU
      for j in range(iAxis):
         Upper *= self.Size[j]
      for j in range(iAxis+1,len(self.Size)):
         Lower *= self.Size[j]
      Shape = (Lower, self.Size[iAxis], Upper)
      return Shape

   def _AssertOpShapeOk(self, Op):
      """assert that an input array has a shape compatible with an periodically
      symmetric operator for this super-cell; that is, either nUnitCells x nSitesU x nSitesU
      or nSitesS x nSitesU."""
      nSitesS = self.nSites
      nSitesU = self.nSitesU
      nGridK = nUnitCells = self.nUnitCells
      assert(nUnitCells * nSitesU == nSitesS)
      assert(Op.shape == (nUnitCells*nSitesU, nSitesU) or
             Op.shape == (nUnitCells,nSitesU,nSitesU))

   def TransformToK(self, Inp):
      """Transform an operator in packed symmetry-unique form (t-space) to
      the symmetry adapted basis (k-space);
      Args:
         - Op: An array of form nUnitCells x nSitesU x nSitesU
               where nSitesU is the number of sites in the unit cell.
               This describes symmetry-unique subset of the operator.
      Returns:
         - An array of form nGridK x nSitesU x nSitesU,
           where nGridK == nUnitCells is the number of k-points in the
           super-cell. Currently, this is always equal to the total number
           of unit-cells, as this allows for one-to-one mappings between
           the two spaces.
      """
      self._AssertOpShapeOk(Inp)

      Out = (self.nUnitCells+0.j)*Inp
      for i in range(len(self.Size)):
         Shape = self._OpShapeForAxis(i)
         Ev = self.kEvs[i]
         Out = (Ev[0,:].reshape(1,self.Size[i],1) *
                np.einsum("tk,ltu->lku", np.conj(Ev), Out.reshape(Shape)))
      return Out.reshape(Inp.shape)

      # hm...
      #if Inp_.dtype in (float,np.float64):
         #assert(np.allclose(Out.imag,0.0))
         #Out = Out.real

      # ^- some notes on what this actually does:
      #
      #  The full transformation is given by
      #
      #     t[K,r,s] = C^H[T,K] t[T,S,r,s] C[S,K]
      #              = C^H[T,K] t[T-S,r,s] C[S,K]
      #
      #  where t are nGridK x nSitesU x nSitesU matrices, R,S are real- space
      #  vectors, and r,s, label sites in the unit cell. The C are the
      #  translation eigenvectors. There is only one K on the lhs, because t is
      #  totally symmetric, and K labels translation irreps.
      #
      #  One must now note that for any given K, the outer product
      #
      #     ProjK[T,S,K] := C^H[T,K] C[S,K]
      #
      #  has the same symmetry as t itself, namely: it depends only on
      #  T-S, not on T and S individually. We can thus in
      #
      #     t[K,r,s] = C^H[T,K] t[T-S,r,s] C[S,K]
      #              = ProjK[T-S,K] t[T-S,r,s]
      #
      #  replace the summation over T or S by a multiplication with nUnitCells,
      #  and have T or S fixed to the 0-vector:
      #
      #     t[K,r,s] = C^H[0,K] t[-S,r,s] C[S,K]
      #              = (C^H[T,K] t[T,r,s]) C[0,K]
      #
      #  We thus only need one actual matrix multiplication over nSitesS
      #  dimension of a sym-uq matrix, and the other transformation is taken
      #  care of by multiplying with a unit-size phase in C[0,K] (which may
      #  be non-1, however)


   def TransformToT(self, Inp):
      """Inverse of self.TransformToK."""
      self._AssertOpShapeOk(Inp)

      Out = Inp
      for i in range(len(self.Size)):
         Shape = self._OpShapeForAxis(i)
         Ev = self.kEvs[i]
         Out = np.einsum("tk,lku->ltu", Ev,
                 np.conj(Ev[0,:].reshape(1,self.Size[i],1)) * Out.reshape(Shape))
      # hm...
      if Inp.dtype in (float,np.float64):
         assert(np.allclose(Out.imag,0.0))
         Out = Out.real
      if Out.dtype in (complex,np.complex128):
         if np.allclose(Out.imag,0.0):
            Out = Out.real
      return Out.reshape(Inp.shape)

   def TransformToEmb(self, OrbL, InpT, OrbR):
      """return Op[k,l] = OrbL[r,k]^* Op[r,s] OrbR[s,l], where
      Op: T x U x U is a periodically symmetric operator expressed
      in the super-cell x unit-cell basis, and OrbL/OrbR are one-particle
      transformations with the super-cell basis as row dimension."""
      nSitesU, nSitesS = self.nSitesU, self.nSites
      assert(InpT.shape == (self.nUnitCells, nSitesU, nSitesU) or
             InpT.shape == (nSitesS, nSitesU))
      assert(OrbL.shape[0] == nSitesS and OrbR.shape[0] == nSitesS)
      Op = InpT.reshape((self.nUnitCells, nSitesU, nSitesU))
      def FindRequiredTs(Op):
         # make a list of Ts we need to consider in the transformation.
         # if Op[iT,:,:] is zero, then we can skip the T.
         iTs = []
         for iT in range(self.nUnitCells):
            if not np.allclose(Op[iT,:,:], 0.):
               iTs.append(iT)
         return np.array(iTs)
      iOpTs = FindRequiredTs(Op)
      #print "iTs with support: %i of %i:  %s" % (len(iOpTs), self.nUnitCells, iOpTs)

      Out = np.zeros((OrbL.shape[1],OrbR.shape[1]),Op.dtype)
      if len(iOpTs) >= .3 * self.nUnitCells:
         Op = InpT.reshape((nSitesS, nSitesU))
         for (i,ijk) in enumerate(self.iTs):
            RowI = self._ConstructTranslatedRows(Op, ijk).T
            # ^- that's the Op[nSitesU*i:nSitesU*(i+1),:] subset of the
            #    operator expressed in the full nSitesS x nSitesS basis.
            #    (see FSuperCell._ExpandFull)
            Out += mdot(np.conj(OrbL[nSitesU*i:nSitesU*(i+1),:].T),
                        RowI, OrbR)
      elif 1:
         Op = Op[iOpTs,:,:]
         OrbR = OrbR.reshape((self.nUnitCells, nSitesU, OrbR.shape[1]))
         for (i,ijk) in enumerate(self.iTs):
            PF,I = self._MakeUnitCellPhasesForT_Restricted(ijk,iOpTs)
            RowIR = np.einsum("t,trs,trl->sl",PF,Op,OrbR[I,:,:])
            Out += np.einsum("rk,rl->kl",np.conj(OrbL[nSitesU*i:nSitesU*(i+1),:]), RowIR)
      else:
         nT = len(iOpTs)
         Op = Op[iOpTs,:,:].reshape((nT*nSitesU,nSitesU))
         nOrbR = OrbR.shape[1]
         OrbR = OrbR.reshape((self.nUnitCells, nSitesU, nOrbR))
         PermShape = (nT*nSitesU,nOrbR)
         PfShape = (nT,1,1)
         for (i,ijk) in enumerate(self.iTs):
            PF,I = self._MakeUnitCellPhasesForT_Restricted(ijk,iOpTs)
            RowIR = np.dot(Op.T, (PF.reshape(PfShape)*OrbR[I,:,:]).reshape(PermShape))
            Out += np.dot(np.conj(OrbL[nSitesU*i:nSitesU*(i+1),:].T), RowIR)
      return Out

   def CalcNonDegenerateOccupations(self, ThrDeg):
      """return a list of all occupation numbers for the current model which
      will result in a non-degenerate ground-state at the tight-binding level"""
      CoreH_T = self.MakeTijMatrix(self.Sites, self.UnitCell)
      CoreH_T = CoreH_T.reshape((self.nUnitCells, self.nSitesU, self.nSitesU))
      CoreH_K = self.TransformToK(CoreH_T)
      AllEw = []
      for ik in xrange(CoreH_K.shape[0]):
         ew,ev = la.eigh(CoreH_K[ik,:,:])
         AllEw.append(ew)
      AllEw = np.array(AllEw).flatten()
      AllEw.sort()
      iOccs = []
      iOcc = 0
      fLast = 1e99
      while iOcc < len(AllEw):
         if ( np.abs(AllEw[iOcc] - fLast) > ThrDeg ):
            iOccs.append(iOcc)
         fLast = AllEw[iOcc]
         iOcc += 1
      #print "Full spectrum: %s" % AllEw
      #print "Occupations with unique fillings: %s" % iOccs
      return iOccs

   def MakeAsy(self, FileName, Colors={}):
      from textwrap import dedent
      from os import system
      Header = dedent("""\
         from plain_pens access *;
         unitsize(20*mm);
         real site_scale = 0.3;
         pen  site_pen = linewidth(1.1) + rgb(0,0,0);
         pen  bond_pen = linewidth(2.9) + rgb(0,0,0);
         pen  bond_pen_warp = linewidth(2.9) + rgb(1.0,0,0);
         void fld(path o, pen paint, pen stroke) { fill(o, paint); draw(o, stroke); }
         defaultpen(font("OT1","cmss","m","n"));
         pair FontShift = - (0,0.04);""")
      Lines = []
      tij = self.Model.MakeTijMatrix(self.Sites, self.Sites)
      if 0:
         tij_sc = self.MakeTijMatrix(self.Sites, self.Sites)
      else:
         tij_sc = tij
      for (iSiteI,SiteI) in enumerate(self.Sites):
         xyi = SiteI[1]
         for (iSiteJ,SiteJ) in enumerate(self.Sites):
            if iSiteI <= iSiteJ or abs(tij_sc[iSiteI, iSiteJ]) < 1e-5:
               continue
            pen = "bond_pen"
            if ( abs(tij[iSiteI, iSiteJ]) < 1e-5 ):
               # sites connect in super-cell T, but not in raw model -> warp
               # hopping due to symmetrization.
               pen = "bond_pen_warp"
            xyj = SiteJ[1]
            lerp = lambda x,y,f: (1.-f)*x + f*y
            xys = lerp(xyi,xyj,0.0)
            xye = lerp(xyi,xyj,1.0)
            Lines.append("draw((%.5f,%.5f)--(%.5f,%.5f), %s);" %
               (xys[0],xys[1],xye[0],xye[1], pen) )
      for (iSite,Site) in enumerate(self.Sites):
         SiteType, SiteXy = Site
         x,y = SiteXy[0], SiteXy[1]
         c = Colors.get(SiteType, (.88,.88,.88))
         Lines.append("fld(circle((%.5f,%.5f), site_scale), rgb(%.2f,%.2f,%.2f), site_pen);" % (x,y,c[0],c[1],c[2]))
         Lines.append("label(scale(1.8)*baseline('%s'),FontShift+(%.5f,%.5f));" % (Site[0],x,y) )
      with open(FileName, "w") as File:
         File.write(Header + "\n" + "\n".join(Lines))
      Cmd = "asy '%s'" % FileName
      print "!%s" % Cmd
      system(Cmd)

def _TestSuperCells():
   from numpy import set_printoptions, nan
   set_printoptions(precision=5,linewidth=10060,suppress=True,threshold=nan)
   if 1:
      L = 12
      nImp = 1
      Hub1d = FHubbardModel1d(nImp, {'t': 1., 'U': 4.})
      sc = FSuperCell(Hub1d, [L/nImp], [-1])
      sc = FSuperCell(Hub1d, [L/nImp], [np.exp(.221j)])
      #sc = FSuperCell(Hub1d, [L/nImp], [1])
   else:
      #Lx,Ly = 320,120
      # ^- FIXME: putting lots of sites there crashes the program. (out of memory)
      #nImpX,nImpY = 2,1
      PhaseShift = None

      FCls = FTestModel2d
      #FCls = FHubbardModel2dSquare

      Lx,Ly=4,3
      #PhaseShift=[-1.,-1.]   # anti-periodic boundary conditions
      PhaseShift=[1.,1.]     # periodic boundary conditions
      PhaseShift=[np.exp(.2j),np.exp(-.72j)] # wtf-like boundary conditions
      nImpX,nImpY = 2,1
      #nImpX,nImpY = 1,1
      nImp = nImpX * nImpY
      Hub2dSq = FCls((nImpX,nImpY), {'t': 1., 'U': 4.})
      sc = FSuperCell(Hub2dSq, [Lx/nImpX,Ly/nImpY], PhaseShift)
   #sc = FSuperCell(Hub1d, [L/nImp,L/nImp], None)
   #sc = FSuperCell(Hub1d, [L/nImp,L/nImp], None)
   CoreH = sc.MakeTijMatrix()
   #print "CoreH via super cells:\n%s" % repr(CoreH)
   print "CoreH spectrum:\n%s" % la.eigh(CoreH)[0]

   #Orb = sc._MakeRealspaceOrb(np.eye(len(sc.Model.UnitCell)), (0,0))
   #print "RealSpaceOrb k == (0,0):\n%s" % Orb
   #print "norm(orb) == ", la.norm(Orb)
   SymBasis = np.hstack([sc._MakeRealspaceOrb(np.eye(len(sc.Model.UnitCell)), ijk) for ijk in sc.iKs])
   print "symmetry adapted basis:\n%s" % SymBasis
   #print "sc.nUnitCells = ", sc.nUnitCells
   #print "sc.nSites = ", sc.nSites
   #print "SymBasis.shape = ", SymBasis.shape
   print "norm( dot(conj(SymBasis.T),SymBasis) - eye(nSites) ): %s" % (la.norm( np.dot(np.conj(SymBasis.T),SymBasis) - np.eye(sc.nSites) ))

   print "CoreH via super cells:\n%s" % CoreH
   print "CoreH in symmetry basis:\n%s" % (mdot(np.conj(SymBasis.T), CoreH, SymBasis))
   #print "CoreH -- 1site transformed:\n%s" % (mdot(CoreH, SymBasis))
   nSitesS = SymBasis.shape[1]
   nSitesU = nImp
   nUnitCells = nSitesS/nSitesU
   CoreH_K = sc.TransformToK(CoreH[:nSitesS,:nSitesU])
   print "And via .TransformToK:\n%s" % CoreH_K
   CoreH_SymUqBt = sc.TransformToT(CoreH_K)
   print "...back to t-space:\n%s" % CoreH_SymUqBt
   print "difference norm to orig: %.2e" % la.norm(CoreH_SymUqBt - CoreH[:nSitesS,:nSitesU])
   raise SystemExit
   #raise Exception("^- that does not look very diagonal (note: works with pbc and apbc, in both directions)")

   nSitesU = len(sc.Model.UnitCell)
   nSitesS = len(sc.Sites)
   #for ijk in sc.iKs.T:
      #Orb = sc._MakeRealspaceOrb(np.eye(len(sc.Model.UnitCell)), ijk)
      #print "ev x evT for ijk = %s:\n%s" % (ijk, np.dot(Orb,np.conj(Orb.T)))
   CoreH_SymUq = CoreH[:nSitesU,:nSitesS]
   print "symmetry unique part of CoreH:\n%s" % CoreH_SymUq
   CoreH_R = np.zeros((nSitesS,nSitesS),complex)
   for (i,ijk) in enumerate(sc.iTs.T):
      CoreH_R[nSitesU*i:nSitesU*(i+1),:] = sc._ConstructTranslatedRows(CoreH_SymUq, ijk)
   print "CoreH orig:\n%s" % CoreH
   print "CoreH reconstructed:\n%s" % CoreH_R
   print "difference norm to orig: %.2e" % la.norm(CoreH_R - CoreH)




if __name__ == "__main__":
   _TestSuperCells()





