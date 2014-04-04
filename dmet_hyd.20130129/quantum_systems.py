# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Gerald Knizia and Garnet K.-L. Chan, 2012;
# (c) Princeton University, 2012

"""This module provides support for some general and non-general quantum
systems described by a Hamiltonian
   H = \sum_{rs} CoreH[rs] c^r c_s + [1/2] \sum_{rstu} Int2e[rstu] c^r c^s c_u c_t

It provides routines for:

   - representing the system (FSystem) and its state (WfDecl),
   - running electronic structure calculations on them (e.g., RunHf for Hartree-Fock)
   - representing the results (FCalcResult)

Also some special code for constructing DMET-style embeddings is provided.
"""
import itertools as it
itt = it
from numpy import *
from numpy.random import random as uniform_random, seed
from scipy.linalg import *
from diis import FDiisContext
from helpers import *
from settings import *
from os import path
from sys import argv
from copy import copy, deepcopy

class FCalcResult(object):
   def __init__(self, Energy, Orbs, OccNumbers, EmbEnergy=None, EmbElec=None, Ews=None):
      self.Energy = Energy
      self.OccNumbers = OccNumbers
      self.Orbs = Orbs # natural orbitals
      if ( EmbEnergy is not None ):
         assert(EmbElec is not None)
         self.EmbEnergy = EmbEnergy
         self.EmbElec = EmbElec
      if ( Ews is not None ):
         self.Ews = Ews # orbital eigenvalues
   def GetRdm(self):
      return mdot(self.Orbs, diag(self.OccNumbers), self.Orbs.T)

class FHfResult(FCalcResult):
   def __init__(self, Energy, Orbs, nOcc, Fock, **pp):
      OccNumbers = zeros(Orbs.shape[0])
      OccNumbers[:nOcc] = ORB_OCC
      FCalcResult.__init__(self,Energy, Orbs, OccNumbers, **pp)
      self.Fock = Fock
   def GetFock(self):
      return self.Fock

class FWfDecl:
   def __init__(self, nElec=None, Ms2=None, nElecA=None, nElecB=None):
      if ( nElecA is not None ):
         assert(nElecB is not None)
         assert(nElec is None and Ms2 is None)
         nElec = nElecA + nElecB
         Ms2 = nElecA - nElecB
      assert(nElec is not None)
      self.nElec = nElec
      self.Ms2 = [Ms2,0][Ms2 is None]
      self.nElecA = (self.nElec + self.Ms2)/2
      self.nElecB = (self.nElec - self.Ms2)/2
      if WF_TYPE == "RHF":
         Ms2 = self.Ms2
         nElec = self.nElec
         assert(Ms2 == 0) # open-shell RHF not supported here.
         assert((nElec + Ms2) % 2 == 0 and Ms2 >= 0 and Ms2 <= nElec)
         self.nOcc = nElec/2
      else:
         self.nOcc = self.nElecA + self.nElecB

def MakeFockOp(Rdm, CoreH, Int2e_Frs, Int2e_4ix=None, Orb=None, ExchFactor=1.):
   """make closed-shell/spin ensemble Fock matrix:
         f[rs] = h[rs] + j[rs] - .5 k[rs].
   where j[rs] = (rs|tu) rdm[tu]
   and   k[rs] = (rs|tu) rdm[su]
   Return (f, j, k)."""
   assert(Int2e_Frs is None or Int2e_4ix is None)
   if Int2e_Frs is not None:
      if 0:
         j = einsum("rs,Frs,Ftu", Rdm, Int2e_Frs, Int2e_Frs)
         k = einsum("rt,Frs,Ftu", Rdm, Int2e_Frs, Int2e_Frs)
         # ^- I wonder if those do something reasonable...
         # ...Update: no.
      else:
         if Orb is None:
            OccNum,Orb = eigh(-Rdm); OccNum *= -1.
         else:
            OccNum = diag(mdot(Orb.T, Rdm, Orb))

         # filter out occupation numbers slightly broken due to rounding.
         OccNum[OccNum < 0.] = 0.
         OccNum[OccNum > ORB_OCC] = ORB_OCC

         nFit, nOrb = Int2e_Frs.shape[0], Int2e_Frs.shape[1]

         # make coulomb.
         gamma = einsum("Frs,rs",Int2e_Frs, Rdm)
         j = einsum("Frs,F",Int2e_Frs,gamma)

         # make exchange -- first need to find the occupied orbitals.
         # (they can be non-standard, e.g., in atomic density initial guess)
         nOcc = len([o for o in OccNum if abs(o-ORB_OCC)<1e-8])
         OccOrb = Orb[:,:nOcc] * OccNum[:nOcc]**.5
         k1 = einsum("Frs,si->Fir", Int2e_Frs, OccOrb).reshape((nFit * nOcc, nOrb))
         k = dot(k1.T, k1)

         Fock = CoreH + j - (ExchFactor/ORB_OCC) * k
         return Fock, j, k
   else:
      assert(Int2e_4ix is not None)
      j = einsum("rstu,tu",Int2e_4ix, Rdm)
      k = einsum("rtsu,tu",Int2e_4ix, Rdm)
   Fock = CoreH + j - (ExchFactor/ORB_OCC) * k
   return Fock, j, k

class ERhfConvergenceFail(Exception):
   def __init__(self, s):
      Exception.__init__(self, s)
   def __str__(self):
      return Exception.__str__(self)

class FSystem(object):
   def __init__(self, CoreEnergy, CoreH, Int2e_Frs, WfDecl, Int2e_4ix=None,
         CoreFockV=None, FockGuess=None, BasisDesc=None):
      """represents a subset of an electronic system (including the full
      system), including its Hamiltonian and wave function declaration.

      Components of H:
         - CoreEnergy: Additive constant to total electronic energy
         - CoreH: 1e matrix elements, nOrb x nOrb matrix
         - 2e integrals can be either given in factorized DF form (normally
           the better way of handling them) or directly as a four-index
           array. The latter form is required if the integral matrix is not
           positive definite. /One/ of the following two should be given:
            o Int2e_Frs: 2e matrix elements in DF/CD form: (F|rs). Actual
              2e matrix elements are then (rs|tu) = \sum_F (rs|F)(F|tu)
              (Mulliken notation).
            o Int2e_4ix: 2e matrix elements in unpacked form:
              nOrb x nOrb x nOrb x nOrb array.
      Other arguments:
         - CoreFockV: Core potential (j/k only) of non-system electrons.
           Not used in WF optimization, but used for calculating
           interaction between special sites and environment.
         - FockGuess: If provided, optionally used for initial guess
           of Fock matrix in RHF or of diagonalization basis in CI.
         - BasisDesc: If given, describes the basis functions of the
           system. May be used to keep track of what is what or to
           form initial guess of Fock matrix (atomic density guess)
      """
      self.CoreH = CoreH
      self.WfDecl = WfDecl
      self.CoreEnergy = CoreEnergy
      assert(Int2e_Frs is None or Int2e_4ix is None)
      self.Int2e_Frs = Int2e_Frs
      self.Int2e_4ix = Int2e_4ix
      self.CoreFockV = CoreFockV

      self.CoreH_Unpatched = 1. * self.CoreH # TODO: improve elegance...
      self.FockGuess = FockGuess

      self.iHfRun = 0

      assert(len(CoreH.shape) == 2 and CoreH.shape[0] == CoreH.shape[1])
      nOrb = CoreH.shape[0]
      self.nOrb = nOrb
      if ( WF_TYPE == "RHF" ):
         self.nSpatialOrb = nOrb
      else:
         assert(WF_TYPE == "UHF" and nOrb%2 == 0)
         self.nSpatialOrb = nOrb/2

      if ( Int2e_Frs is not None ):
         assert(len(Int2e_Frs.shape) == 3 and Int2e_Frs.shape[1] == Int2e_Frs.shape[2] == nOrb)
         self.nFit = Int2e_Frs.shape[0]
      else:
         assert(Int2e_4ix.shape == (nOrb,nOrb,nOrb,nOrb))
         self.nFit = None
      self.BasisDesc = BasisDesc
      pass

   def MakeImpPrjSystem(self, SpecialSites):
      """construct a new system object which represents the fake
      Hamiltonian one can use in <x|FakeH|y> to calculate the energy
      within the special sites, and their interaction energy with
      the rest of the system.

      Note that this system will generally need to have a real 4ix 2e
      interaction, because the projected Int2e is not positive
      semidefinite in general."""
      OtherSites = array([i for i in xrange(self.nOrb) if i not in SpecialSites], int)

      # Note on the the 2e integrals:
      # The following two variants give the same result. The first is
      # the reference version: The energy formula for the full system is
      #
      #     E = h[rs] rdm1[rs] + .5 * (V[rstu] rdm2[rstu])
      #
      # In order to calulate the energy associated with the impurity, we
      # /fix/ r to r\in imp. This will return only the energy within the
      # impurity itself and for the impurity/rest interaction (one
      # time). Thus we get:
      #
      #     E_{imp,imp-bath} = h[0s] rdm1[0s] + .5 * (V[0stu] rdm2[0stu])
      #
      # however, since both V and the 2rdm have the same 8-fold
      # symmetry, we can get the same result by setting the non- present
      # rows of V to zero and symmetrizing the resulting V. This is what
      # is happening here. The second version shows that this results in
      # various factors of 0.,.25,.5,.75,and 1. depending on the number
      # of sites in the impurity, in the resulting 2e integrals.
      nOrb = self.nOrb
      V = 1.*self.MakeOrGet_Int2e_4ix()

      if 1:
         V[OtherSites,:,:,:] = 0
         V = 1/8.*(V.transpose([0,1,2,3]) + V.transpose([2,3,0,1]) +\
                   V.transpose([1,0,2,3]) + V.transpose([3,2,0,1]) +\
                   V.transpose([0,1,3,2]) + V.transpose([2,3,1,0]) +\
                   V.transpose([1,0,3,2]) + V.transpose([3,2,1,0]))
      else:
         f = array([0., .25, .5, .75, 1.])
         for i,j,k,l in it.product(xrange(nOrb),xrange(nOrb),xrange(nOrb),xrange(nOrb)):
            nsys = int(i in SpecialSites) + int(j in SpecialSites) +\
                   int(k in SpecialSites) + int(l in SpecialSites)
            V[i,j,k,l] *= f[nsys]

      # Note on the core Fock contribution:
      # If applied to the full system (not the embedded system), the 2e
      # argument can also be used to derive the scaled core Fock
      # contribution. Note that when involving environment electrons
      # (i.e., neither sys nor bath), the 2-RDM factorizes into an
      # impurity 1-RDM and an environment determinant 2-RDM:
      #    rdm2[i,r,s,j] = <c^i c^r c_s c_j>_full
      #                  = <c^i c_j>_env <c^r c_s>_imp
      #                  = delta[i,j] <c^r c_s>_imp
      # (with i,j \in env and r,s \in imp, and similar for other cases
      # and the spatial version). Thus, in the actual interacting 2e
      # integrals there are only contributions (rs|ii) and (ri|si) with
      # i \in env and r,s \in imp, and both go into the 1e-part of the
      # embedded system. Since there are never three indices in the
      # impurity, the only scaling factors occuring are .5,.25, and 0.,
      # as we have above.

      #... hmhm.. is that -^ right? Are there really no (virtual?) environment
          #contributions with factor 0.75?

      SubCoreH = 1. * self.CoreH
      if ( self.CoreFockV is not None ):
         # j/k interaction with the non-system electrons is also
         # an interaction term, and must be scaled by 1/2!
         # currently it's fully included in CoreH.
         SubCoreH -= .5 * self.CoreFockV

      for (i,j) in it.product(OtherSites, OtherSites):
         SubCoreH[i,j] = 0.
      for (i,j) in it.product(SpecialSites, OtherSites):
         SubCoreH[i,j] *= 0.5
         SubCoreH[j,i] *= 0.5

      return FSystem(self.CoreEnergy, SubCoreH, None, self.WfDecl, Int2e_4ix = V)

   def Run(self, Method, Print, SpecialSites = None):
      """If SpecialSites is given, also calculate the embedded energy of the
      special sites, which consists of the energy within the special sites
      themselves and their interaction energy with the rest."""
      if Method == "RHF":
         return self.RunHf(Print, SpecialSites)
      elif Method == "FCI" or Method.startswith('CI(') or Method.startswith('CC('):
         return self.RunFci0(Print, SpecialSites, CiMethod = Method)
      elif Method == "MP2" or Method.startswith('CI(') or Method.startswith('CC('):
         return self.RunMp2(Print, SpecialSites)
      else:
         raise Exception("Electronic structure method not recognized: '%s'" % Method)

   def CalcHfEnergy(self, Rdm):
      # calculate the RHF energy for a given density matrix.
      Fock, j, k = MakeFockOp(Rdm, self.CoreH, self.Int2e_Frs, self.Int2e_4ix)
      Energy = 0.5 * (dot2(Fock, Rdm) + dot2(self.CoreH, Rdm)) + self.CoreEnergy
      return Energy

   def ProjectToIrrep0(self, M):
      assert(hasattr(self.WfDecl,"SoBasis"))
      SoBasis = self.WfDecl.SoBasis
      Out = zeros_like(M)
      for IrpBasis in SoBasis:
         Out += mdot(IrpBasis, mdot(IrpBasis.T, M, IrpBasis), IrpBasis.T)
      return Out

   def RunHf(self, Print, SpecialSites = None):
      self.iHfRun += 1
      p = Print
      nOrb = self.nOrb
      if ( self.FockGuess is None ):
         if ( self.BasisDesc is None ):
            # coreH guess
            Fock = self.CoreH
         else:
            # atomic density guess
            DiagonalOcc = [self.BasisDesc[i].AtomicOccupancy for i in range(nOrb)]
            Rdm = diag(DiagonalOcc)
            Fock, j, k = MakeFockOp(Rdm, self.CoreH, self.Int2e_Frs,
               Int2e_4ix=self.Int2e_4ix, ExchFactor=0.5)
            pass
      else:
         Fock = self.FockGuess

      #Fock = 1.*self.CoreH_Unpatched
      if hasattr(self.WfDecl,"SoBasis"): Fock = self.ProjectToIrrep0(Fock)
      if hasattr(self.WfDecl,"SoBasis"): self.CoreH = self.ProjectToIrrep0(self.CoreH)

      LastEnergy = 0.
      nOcc = self.WfDecl.nOcc
      c0 = 0.
      if p:
         print "*** %sHARTREE-FOCK\n" % (["","UNRESTRICTED "][WF_TYPE=="UHF"])
         print " %-35s%5i FUNCTIONS" % ("ORBITAL BASIS:", nOrb)
         print " %-35s%5i FUNCTIONS" % ("FITTING BASIS:", self.Int2e_Frs.shape[0])
         print

      def UpdateOrbitals(Fock,RefRdm,iSpin=None):
         Ew,Orbs = eigh(Fock)
         if RefRdm is not None:
            # freeze occupation pattern. Occupy orbitals with maximum
            # overlap with previous occupied space.
            PrevOcc = diag(mdot(Orbs.T, RefRdm, Orbs))
            #if iSpin==0: print "\noverlap with previous orbtials:\nA%s" % PrevOcc
            #if iSpin==1: print "B%s\n" % PrevOcc
            iOcc = PrevOcc.argsort()[::-1]
            if 1:
               Ew = Ew[iOcc]
               Orbs = Orbs[:,iOcc]
         return Ew,Orbs

      Converged = False
      RefRdm = None
      dc = FDiisContext(12)
      if p: print "   {:^5s} {:^14s} {:^14s} {:^11s} {:>4s}".format(
         "ITER.","ENERGY","ENERGY CHANGE", "GRAD", "DIIS")
      LastRdm = None
      for it in range(MAXIT_HF):
         # make new orbitals and new rdm.
         if WF_TYPE == "RHF":
            Ew,Orbs = UpdateOrbitals(Fock,RefRdm)
            #PrintMatrix("Fock spectrum", Ew)
         else:
            ESC = ExtractSpinComp
            EwA,OrbsA = UpdateOrbitals(ESC(Fock,0), ESC(RefRdm,0), 0)
            EwB,OrbsB = UpdateOrbitals(ESC(Fock,1), ESC(RefRdm,1), 1)
            Ew = CombineSpinComps(EwA,EwB)
            Orbs = zeros((nOrb,nOrb))
            Orbs[ ::2, ::2] = OrbsA[:,:]
            Orbs[1::2,1::2] = OrbsB[:,:]

         Rdm = ORB_OCC * dot(Orbs[:,:nOcc], Orbs[:,:nOcc].T)
         if hasattr(self.WfDecl,"SoBasis"): Rdm = self.ProjectToIrrep0(Rdm)

         # make closed-shell Fock matrix:
         #   f[rs] = h[rs] + j[rs] - .5 k[rs]
         Fock, j, k = MakeFockOp(Rdm, self.CoreH, self.Int2e_Frs,
                                 Int2e_4ix=self.Int2e_4ix, Orb=Orbs)
         Energy = 0.5 * (dot2(Fock, Rdm) + dot2(self.CoreH, Rdm)) + self.CoreEnergy
         if hasattr(self.WfDecl,"SoBasis"): Fock = self.ProjectToIrrep0(Fock)

         # calculate orbital gradient and energy
         OrbGrad = dot(Fock,Rdm)
         if WF_TYPE == "UHF":
            OrbGrad[ ::2,1::2] = 0
            OrbGrad[1::2,0::2] = 0
         OrbGrad = OrbGrad - OrbGrad.T
         #print "Orbital gradient:\n%s"%OrbGrad
         # ^- fun fact: "OrbGrad -= OrbGrad.T" is an unsafe operation
         #    which, depending on circumstances, might or might not produce
         #    the expected result (well.. kinda obvious in hindsight).
         if 0:
            PrintMatrix("Orbital energies:", eigh(Fock)[0])
            Fai = mdot(Orbs[:,nOcc:].T, Fock, Orbs[:,:nOcc])
            PrintMatrix("F[a,i] Fock matrix elements", Fai)
            OrbGrad = mdot(Orbs[:,nOcc:], Fai, Orbs[:,:nOcc].T)

         fOrbGrad = norm(OrbGrad.flatten())
         dEnergy = Energy - LastEnergy
         LastEnergy = Energy

         if p: print " {:5d} {:14.8f} {:+14.8f} {:11.2e} {:>6s}".format(
            it+1, Energy, dEnergy, fOrbGrad, dc)

         if (fOrbGrad < THRORB and abs(dEnergy) < THRDEN):
            Converged = True
            break
         if (it == MAXIT_HF - 1):
            break # not converged.

         # extrapolate
         SkipDiis = False
         Fock, OrbGrad, c0 = dc.Apply(Fock, OrbGrad, Skip=SkipDiis)
      #if ( self.iHfRun == 2 ): raise SystemExit

      if (not Converged and MAXIT_HF != 1):
         if 1:
            print "CoreH:\n%s" % self.CoreH
            print "Vloc:\n%s" % (self.CoreH - self.CoreH_Unpatched)
            print "Fock:\n%s" % Fock
            print "Spectrum:\n %s" % eigh(Fock)[0]
            print "RDM:\n%s" % Rdm
         ErrMsg = "%s failed to converge."\
                  "  NIT =%4i  GRAD=%8.2e  DEN=%8.2e" % (WF_TYPE, it+1, fOrbGrad, dEnergy)
         print "WARNING: %s" % ErrMsg
         if 1:
            print "Failed FSystem object stored to file 'brrrk' for investigation."
            import pickle
            pickle.dump(self,open('brrrk','w'))
         raise ERhfConvergenceFail(ErrMsg)
      elif p:
         print
         print pResultFmt % ("1e energy", dot2(Rdm, self.CoreH))
         print pResultFmt % ("2e energy", .5*dot2(Rdm, Fock - self.CoreH))
         if (self.CoreEnergy != 0.0):
            print
            print pResultFmt % ("Core energy", self.CoreEnergy)
            print pResultFmt % ("Electronic energy", Energy - self.CoreEnergy)
         print pResultFmt % ("Hartree-Fock energy", Energy)
         print

      EmbEnergy = None
      EmbElec = None
      if (SpecialSites is not None):
         # calculate number of electrons and energy in the embedded system.
         EmbElec = trace(Rdm[:,SpecialSites][SpecialSites,:])
         FakeSys = self.MakeImpPrjSystem(SpecialSites)
         EmbEnergy = FakeSys.CalcHfEnergy(Rdm)
         if p:
            print pResultFmt % ("Embedded subsystem energy", EmbEnergy)
            print pResultFmt % ("Embedded subsystem electrons", EmbElec)
            print
      return FHfResult(Energy, Orbs, nOcc, Fock, EmbElec=EmbElec, EmbEnergy=EmbEnergy, Ews=Ew)

   def RunMp2(self, Print, SpecialSites):
      p = Print
      HfResult = self.RunHf(Print,SpecialSites)

      Orbs = HfResult.Orbs
      Ew = diag(mdot(Orbs.T,HfResult.GetFock(),Orbs))
      nOcc = self.WfDecl.nOcc
      nOrb = self.nOrb
      def MakeMp2T():
         OrbOcc = Orbs[:,:nOcc]; EwOcc = Ew[:nOcc]
         OrbVir = Orbs[:,nOcc:]; EwVir = Ew[nOcc:]
         K_Fri = einsum("Frs,si->Fri", self.Int2e_Frs, OrbOcc)
         K_Fai = einsum("Fri,ra->Fai", K_Fri, OrbVir)
         Kabij = einsum("Fai,Fbj->abij", K_Fai, K_Fai)
         o, v = OrbOcc.shape[1], OrbVir.shape[1]
         Denoms = EwVir.reshape(v,1,1,1) + EwVir.reshape(1,v,1,1) - EwOcc.reshape(1,1,o,1) - EwOcc.reshape(1,1,1,o)
         Tabij = -Kabij / Denoms
         TabijB = 2.*Tabij - Tabij.transpose(1,0,2,3)

         RdmMo = zeros((nOrb,nOrb))
         RdmMo[:nOcc,:nOcc] = 2.*eye(nOcc) - 2.*einsum("abmi,abni",TabijB,Tabij)
         RdmMo[nOcc:,nOcc:] = 2.*einsum("agij,ahij",Tabij,TabijB)
         Rdm = mdot(Orbs, RdmMo, Orbs.T)

         Mp2Energy = dot2(TabijB,Kabij)
         return einsum("abij,ra,sb,ti,uj->rstu",TabijB,OrbVir,OrbVir,OrbOcc,OrbOcc).transpose(0,2,1,3), Rdm, Mp2Energy
      V = self.MakeOrGet_Int2e_4ix()
      Trstu, Rdm, Mp2Energy = MakeMp2T()

      if (SpecialSites is not None):
         # calculate number of electrons and energy in the embedded system.
         EmbElec = trace(Rdm[:,SpecialSites][SpecialSites,:])
         FakeSys = self.MakeImpPrjSystem(SpecialSites)
         #EmbEnergy = FakeSys.CoreEnergy + dot2(FakeSys.CoreH, Rdm) + dot2(FakeSys.Int2e_4ix, Trstu)
         #print pResultFmt % ("2*Embedded subsystem 1e energy (M)", 2*dot2(Rdm, FakeSys.CoreH))
         EmbEnergy = HfResult.EmbEnergy + dot2(FakeSys.Int2e_4ix, Trstu)
         #EmbEnergy += dot2(FakeSys.CoreH, Rdm - HfResult.GetRdm())
         # ^- should this be here? This form of MP2 is really not supposed to have singles contributions.
         if p:
            print pResultFmt % ("Embedded subsystem energy", EmbEnergy)
            print pResultFmt % ("Embedded subsystem electrons", EmbElec)
            print

      OccNumbers,Orbs = eigh(-Rdm); OccNumbers *= -1.
      if p: print " %-32s   %-s" % ("Natural occupation", " ".join(["%7.4f" % o for o in OccNumbers]))
      return FCalcResult(HfResult.Energy + Mp2Energy, Orbs, OccNumbers, EmbElec=EmbElec, EmbEnergy=EmbEnergy)


   def MakeOrGet_Int2e_4ix(self):
      if ( self.Int2e_4ix is not None ):
         return self.Int2e_4ix
      else:
         assert(self.Int2e_Frs is not None)
         return einsum("Frs,Fut",self.Int2e_Frs, self.Int2e_Frs)

   def RunFci(self, Print, SpecialSites = None):
      return self.RunFci0(Print, SpecialSites)

   def RunFci0(self, Print, SpecialSites = None, CiMethod="FCI"):
      import fci_iface as fi
      from os import remove
      #print "via FCI: ",;  DummyHfResult = self.RunHf(Print,SpecialSites) # FIXME: remove this
      p = Print
      def DeleteFciVectorFile():
         # remove left-over FCI vector if present
         if fi.FileNameFciVec:
            try: remove(fi.FileNameFciVec)
            except OSError,e: pass
      DeleteFciVectorFile()

      nSysTrace=None
      if ( SpecialSites is not None ):
         nSysTrace = len(SpecialSites)
         assert((array(SpecialSites) == arange(0, len(SpecialSites))).all())
      HL = fi.FindFciSolution(self.CoreH, U=1., nSys_=self.nOrb, # <- both unused.
         nElec=self.WfDecl.nElec, CiMethod=CiMethod,Int2e_=self.MakeOrGet_Int2e_4ix(),
         nSysTrace=nSysTrace, IntType=WF_TYPE, MakeRdm2=False, DiagBasis=eigh(self.FockGuess)[1])
      if ( IPRINT >= 4 ):
         print HL["Output"] + "\n"
      Energy, rdm = HL["Energy"], HL["Rdm1"]
      OccNumbers,Orbs = eigh(-rdm); OccNumbers *= -1.
      EmbEnergy = None
      EmbElec = None
      if p:
         print " %-32s   %-s" % ("Natural occupation", " ".join(["%7.4f" % o for o in OccNumbers]))
         print pResultFmt % ("Full-CI energy", Energy)
         print
      if ( SpecialSites is not None ):
         EmbElec = trace(rdm[:,SpecialSites][SpecialSites,:])
         FakeSys = self.MakeImpPrjSystem(SpecialSites)
         if 0:
            # this should just read the fci vector from file and calculate its energy.
            # ...hmpf. need to make sure that both use the same diagonalization basis.
            HL1 = fi.FindFciSolution(FakeSys.CoreH, 1., FakeSys.nOrb/2,
               FakeSys.WfDecl.nElec, "FCI",Int2e_=FakeSys.MakeOrGet_Int2e_4ix(),MaxIt=1)
            EmbEnergy = HL1["Energy"]
         else:
            # should also do the trick... 2e energy calculated by FCI program itself.
            EmbEnergy = dot2(rdm, FakeSys.CoreH) + HL["EpTrace"] + FakeSys.CoreEnergy

         if p:
            print pResultFmt % ("Embedded subsystem energy", EmbEnergy)
            print pResultFmt % ("Embedded subsystem electrons", EmbElec)
            print

      if 0:
         print HL["Output"]
         print HL["Cmd"]
         print "waiting for [enter] or [Ctrl+C]..."
         raw_input()
      DeleteFciVectorFile()
      return FCalcResult(Energy, Orbs, OccNumbers, EmbElec=EmbElec, EmbEnergy=EmbEnergy)

   def CalcFciEnergy1(self, CiVector):
      """calculate the FCI energy for a given CI vector."""
      import StringDet as sd
      V = self.MakeOrGet_Int2e_4ix()
      Configs = sd.MakeFullConfigList(self.WfDecl.nElecA, self.WfDecl.nElecB, self.nSpatialOrb)
      assert(len(Configs) == len(CiVector))
      H = sd.MakeDeterminantH(Configs, self.CoreH, V, IntType=WF_TYPE)
      #print "CiVector:\n ", CiVector
      #print "DeterminatH:\n%s"%H
      #print "Spectrum:\n %s"%eigh(H)[0]
      #raise SystemExit
      return einsum("I,IJ,J", CiVector, H, CiVector)

   def RunFci1(self, Print, SpecialSites = None):
      """run fci via python StringDet.py"""
      import StringDet as sd
      p = Print
      # make 2e integrals in unpacked form.
      V = self.MakeOrGet_Int2e_4ix()
      Configs = sd.MakeFullConfigList(self.WfDecl.nElecA, self.WfDecl.nElecB, self.nSpatialOrb)
      if p:
         print "*** FULL-CI\n"
         print " %-35s%5i FUNCTIONS" % ("ORBITAL BASIS:", self.nOrb)
         print " %-35s%5i FUNCTIONS" % ("DETERMINANT BASIS:", len(Configs))
         print

      # form and diagonalize the full Hamiltonian in determinant basis.
      H = sd.MakeDeterminantH(Configs, self.CoreH, V, IntType=WF_TYPE)
      iRoot = 0
      ew,ev = eigh(H)
      Energy = ew[iRoot] + self.CoreEnergy
      CiVector = ev[:,iRoot] # ground state CI vector.
      rdm = sd.Make1Rdm(Configs, CiVector, CiVector, IntType=WF_TYPE)
      OccNumbers,Orbs = eigh(-rdm); OccNumbers *= -1.

      EmbEnergy = None
      EmbElec = None
      if p:
         print " %-32s   %-s" % ("Natural occupation", " ".join(["%7.4f" % o for o in OccNumbers]))
         print pResultFmt % ("Full-CI energy", Energy)
         print
      if ( SpecialSites is not None ):
         EmbElec = trace(rdm[:,SpecialSites][SpecialSites,:])
         FakeSys = self.MakeImpPrjSystem(SpecialSites)
         EmbEnergy = FakeSys.CalcFciEnergy1(CiVector)
         if p:
            print pResultFmt % ("Embedded subsystem energy", EmbEnergy)
            print pResultFmt % ("Embedded subsystem electrons", EmbElec)
            print
      return FCalcResult(Energy, Orbs, OccNumbers, EmbElec=EmbElec, EmbEnergy=EmbEnergy)

def ReRunFailedHf():
   import pickle
   sys = pickle.load(open('brrrk','r'))
   sys.RunHf(Print=True)
   raise SystemExit
if __name__ == "__main__":
   ReRunFailedHf()




