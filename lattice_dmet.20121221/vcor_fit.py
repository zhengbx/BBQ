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

"""provides support for calculating potential matrix updates based on
a high-level density matrix."""
from helpers import mdot
import numpy as np
from numpy import dot, zeros, zeros_like, array, diag, einsum
import scipy.linalg as la

def FitVcorComponent(EmbFock, nImp, RdmHl, VLOC_TYPE, VLOC_FIT_TYPE):
   """Fit an update matrix for a given RdmHl. This routine does not
   know about spin: It assumes the all orbitals are SINGLY OCCUPIED.
   (that means, in RHF mode, you must supply RdmHl/2 are argument,
    while in UHF mode, you supply matrices for the individual spin
    components)
   Inputs:
   - EmbFock: current Fock matrix
   - nImp: Number of impurity sites on which to put the potential.
           Determines output dimension; output is nImp x nImp matrix.
   - RdmHl: high-level embedded density matrix. Must have same shape
            as EmbFock.
   Returns nImp x nImp matrix dVcor.
   """
   assert(EmbFock.shape[0] == EmbFock.shape[1])
   assert(EmbFock.shape == RdmHl.shape)
   if VLOC_TYPE == "Diagonal":
      VlocSize = nImp
   else:
      assert(VLOC_TYPE == "Local")
      VlocSize = (nImp*(nImp+1))/2

   # nEmb: total number of orbitals in embedded system.
   nEmb = EmbFock.shape[0]

   def MakePotentialPatchMatrix(VlocSize):
      # make a Matrix nEmb x nEmb x VlocSize, mapping vloc
      # components (last dimension) to the potential matrix they
      # produce (first two dimensions)
      h0_for_vloc = zeros((nEmb,nEmb,VlocSize))
      ij = 0
      for i in range(nImp):
         for j in range(i+1):
            if ( VLOC_TYPE == "Diagonal" and i != j ):
               continue
            h0_for_vloc[i,j,ij] = 1.
            h0_for_vloc[j,i,ij] = 1.
            ij += 1
      return h0_for_vloc.reshape((nEmb*nEmb, VlocSize))
   h0_for_vloc = MakePotentialPatchMatrix(VlocSize)

   def GetImpRdmFitAnalytical(RdmHl, nImp, VLOC_TYPE, VLOC_FIT_TYPE):
      # this does the same as resmin, but
      # solves the target fitting equation analyically.
      # note that we can setup an equivalent routine for fullrdm fit.
      assert(VLOC_TYPE == "Local")
      # ^- we could do Diagonal/ImpRdm by making an idempotent
      # approximation to a RDM obtained by taking the previous HF
      # RDM and just fixing the diagonal entries to the CI RDM.

      def MakeIdempotentApproxRdm(TargetRdm):
         # make a idempotent density matrix which exactly
         # coincides with TargetRdm on the imp x imp subset.
         A = TargetRdm[:nImp,:nImp]
         C = TargetRdm[nImp:,nImp:]
         assert(A.shape == C.shape)
         C = np.eye(A.shape[0]) - A
         # ^-
         # WARNING: at this moment this does assume that the embedding
         #   basis is constructed as full Smh basis, i.e., that the hole
         #   state corresponding to imp site i is at nImp + i. For other
         #   cases we need additional rotations on C (not simply overwrite
         #   it) and fix the eigenvalues in some other way.
         def MakeRho(A,B,C):
            return np.vstack([np.hstack([A,B]),np.hstack([np.conj(B.T),C])])

         L = nImp
         ewA,evA = la.eigh(A)
         ewC,evC = la.eigh(C)
         basis = MakeRho(evA,zeros_like(A),evC)
         R = mdot(basis.T, TargetRdm, basis)
         M = 0.*basis
         M[:L,:L] = diag(ewA)
         M[L:,L:] = diag(ewC)
         ewB = (ewA - ewA**2)**.5
         # ^- we can still put in random phases here and retain an
         #    idempotent target matrix. we get them from the
         #    transformation of the target density matrix, hopefully
         #    in such a way that also the off-diagonals are not way off,
         #    even if we are only fitting on ImpRdm instead of FullRdm.
         #    (see idempotency_gen.py)
         #
         #    This freedom may be related to the freedom in choosing orbital
         #    occupation.
         ewB_ref = np.diag(R[nImp:,:nImp][::-1,:])
         ewB *= np.sign(ewB_ref)
         M[L:,:L] = np.diag(ewB)[::-1,:]
         M[:L,L:] = np.conj(M[L:,:L]).T
         rho = dot(basis,np.dot(M,basis.T))
         assert(np.allclose(np.dot(rho,rho)-rho,0.0))
         # ^- will work for HF in HF, but not for FCI in HF
         return rho

      def MakeVlocForRdm(NewRdm, nOcc, EmbFock, h0_for_vloc):
         # okay... we have the target eigensystem. We now still
         # need an operator EmbFock + vloc which generates it as
         # a ground state density. The condition here is that
         #
         #    dot(h0_for_vloc, vloc1)  =: K v
         #
         # has the same eigensystem as the one we're calculating,
         # or, equivalently, that it commutes with the rdm we have.
         # We assert that the resulting Fock matrix does not have
         # components mixing occupied and virtual orbitals, which is
         # equivalent.
         # Note that we cannot assume that F is diagonal in the target rdm
         # eigenvectors, because the eigensystem is degenerate
         # (has only 0 and 1 evs)
         #
         # Also note that in general vloc has not enough variational
         # freedom to represent arbitrary eigensystems, so technically we
         # are returning the best approximation which is possible.
         ew,ev = la.eigh(NewRdm)
         O,V = ev[:,nOcc:], ev[:,:nOcc]
         lhs = mdot(O.T,EmbFock,V).flatten()
         VlocSize = h0_for_vloc.shape[-1]
         h0_for_vloc_reshape = h0_for_vloc.reshape((nEmb,nEmb,VlocSize))
         K = einsum("ro,rsp,sv->ovp",O,h0_for_vloc_reshape,V).reshape(O.shape[1]*V.shape[1],VlocSize)
         rhs = -la.lstsq(K,lhs)[0]
         return dot(h0_for_vloc_reshape,rhs)

      # get number of orbitals we need to occupy in embedded HF.
      fOcc = np.trace(RdmHl)
      nOcc = int(fOcc+.5)
      assert(abs(fOcc - nOcc) < 1e-8)

      if VLOC_FIT_TYPE == "ImpRdm":
         # make an idempotent density matrix which exactly
         # coincides with RdmHl on the Imp x Imp part.
         FitRdmHf = MakeIdempotentApproxRdm(RdmHl)
      elif VLOC_FIT_TYPE == "FullRdm":
         # make the best idempotent approximation to the FCI
         # RDM. That's simply the matrix with the same natural
         # orbitals, but all occupation numbers rounded to 0
         # (if < .5) or 1 (if > .5).
         if 0:
            OccNum,NatOrb = la.eigh(-RdmHl); OccNum *= -1.
            OccNumHf = zeros_like(OccNum)
            OccNumHf[:nOcc] = 1.
            FitRdmHf = mdot(NatOrb, diag(OccNumHf), NatOrb.T)

         # ^- okay.. that is not *entirely* right:
         #  here we first approximate the target RDM, and then
         #  approximate the potential changes required to most
         #  closely obtain the target RDM as a HF solution.
         FitRdmHf = RdmHl
      else:
         # not supported.
         assert(0)
      def TestMakeVlocForRdm():
         vloc = MakeVlocForRdm(FitRdmHf, nOcc, EmbFock, h0_for_vloc)
         NewF = EmbFock + vloc
         ew,ev = la.eigh(NewF)
         TestRdm = mdot(ev[:,:nOcc], ev[:,:nOcc].T)
         print TestRdm-FitRdmHf
      return MakeVlocForRdm(FitRdmHf, nOcc, EmbFock, h0_for_vloc)

   if VLOC_TYPE == "Local" and VLOC_FIT_TYPE in ["ImpRdm", "FullRdm"]:
      # hm.. this is actually not quite 100% the same as we did before
      # in the FullRdm case. It gives, however, almost the same numbers
      # and is more stable and much more scalable...
      return GetImpRdmFitAnalytical(RdmHl, nImp, VLOC_TYPE, VLOC_FIT_TYPE)
   else:
      from helpers import resmin

      # find number of electrons/orbitals to occupy.
      fElec = np.trace(RdmHl)
      nElec = int(fElec+.5)
      assert(abs(fElec - nElec) < 1e-8)
      nOcc = nElec

      def MakeF(vloc):
         if vloc is None: return EmbFock
         return EmbFock + dot(h0_for_vloc, vloc).reshape((nEmb,nEmb))

      def MakeRdm(vloc):
         F = MakeF(vloc)
         ew,ev = la.eigh(F)
         OrbCl = ev[:,:nOcc]
         return dot(OrbCl, OrbCl.T)

      RdmHf = MakeRdm(vloc=None)

      Occ,NatOrb = la.eigh(-RdmHf); Occ *= -1.
      for (io,o) in enumerate(Occ):
         if (o < 0.): Occ[io] = 0.
      OccOrb = dot(NatOrb, diag(Occ**.5))
      def Err1v(vloc):
         dRdm = (MakeRdm(vloc) - RdmHl)
         if VLOC_TYPE == "Diagonal":
            dRdm = diag(dRdm)
         if VLOC_FIT_TYPE == "FullRdm":
            return dRdm.flatten()
         elif VLOC_FIT_TYPE == "ImpRdm":
            if VLOC_TYPE == "Diagonal":
               return dRdm[:nImp]
            return ToTriangle(dRdm[:nImp,:nImp])
         elif VLOC_FIT_TYPE == "EnvRdm":
            if VLOC_TYPE == "Diagonal":
               return dRdm[nImp:]
            return ToTriangle(dRdm[nImp:,nImp:])
         elif VLOC_FIT_TYPE == "EnvAndCoupling":
            dRdm[:nImp,:nImp] = 0
            return dRdm.flatten()
         else:
            # fit type not recognized.
            assert(0)

      vloc1 = resmin(Err1v, zeros(VlocSize), Algo="Hard")
      vloc = dot(h0_for_vloc, vloc1).reshape((nEmb,nEmb))
      return vloc
