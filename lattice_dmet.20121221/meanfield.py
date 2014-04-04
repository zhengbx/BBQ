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

"""This file contains functions and classes for dealing with the mean-field
of the full lattice[*]. It is not concerned with the mean-field in the embedded
quantum systems.

[*] by that I mean the cyclic super-cell model of the entire lattice.
"""
import numpy as np
from scipy import linalg as la
import itertools as it
from lattice_model import *
from supercell import *
from output import *
from helpers import mdot, dot2, ExtractSpinComp, CombineSpinComps
from diis import FDiisContext
from quantum_systems import FWfDecl

#class FLatticeMfResult(object):
   #def __init__(self, Energy, Orbs, OccNumbers, Fock, Ews, P):
      #assert(Fock.shape == Orbs.shape)
      ## Fock: T x U x U matrix (T: #UnitCells in super-cell, U: nSites in unit-cell)
      #self.Fock = Fock
      ## Orbs: K x U x U matrix (K: number of K-points (==#Unitcells))
      ## Note: orbs are stored in symmetry adapted basis (k-space basis), not
      ##       real-space basis (t-space)
      #self.Orbs = Orbs
      ## OccNumbers: K x U matrix;
      ## [k,r]: occupation number of orbital r in irrep k.
      #self.OccNumbers = OccNumbers
      #self.Energy = Energy
      ## Ews: orbital eigenvalues (must be equal to <r|F|s> where F is the Fock operator)
      ## Same format as OccNumbers.
      #self.Ews = Ews
      #pass

   #def GetRdm(self):
      ##return 2.0 * dot(self.Orbs[:,:nOcc], self.Orbs[:,:nOcc].T)
      ##OccOrbs = 0
      #return mdot(self.Orbs, diag(self.OccNumbers), self.Orbs.T)

class FMeanFieldParams:
   def __init__(self, ITERATIVE="auto", THRORB=1e-8, THRDEN=1e-8,
         THRDEG=1e-5, MAXIT=1024, DIIS_DIM=12, DIIS_START=0, DIIS_THR=1e99):
      self.ThrDen = THRDEN
      self.ThrOrb = THRORB
      self.ThrDeg = THRDEG
      self.MaxIt = MAXIT
      self.DiisDim = DIIS_DIM
      self.DiisStart = DIIS_START
      self.DiisThr = DIIS_THR
      self.Iterative = ITERATIVE


class EHfConvergenceFail(Exception):
   def __init__(self, s):
      Exception.__init__(self, s)
   def __str__(self):
      return Exception.__str__(self)

def ExtractSpinCompTK(M, iComp):
   assert(0<= iComp <= 1)
   if len(M.shape) == 3:
      return M[:, iComp::2,iComp::2]
   elif len(M.shape) == 2:
      return M[:, iComp::2]
   else:
      assert(0)

def CombineSpinCompsTK(M_Alpha, M_Beta):
   assert(M_Alpha.shape == M_Beta.shape)
   Shape = [M_Alpha.shape[0]] + [2*o for o in M_Alpha.shape[1:]]
   M = np.zeros(Shape, M_Alpha.dtype)
   if len(M.shape) == 3:
      M[:, ::2, ::2] = M_Alpha
      M[:,1::2,1::2] = M_Beta
      return M
   elif len(M.shape) == 2:
      M[:, ::2] = M_Alpha
      M[:,1::2] = M_Beta
   else:
      assert(0)
   return M

#def MakeRdm(Ew,Orbs,nElec,OccAvg):
   #from scipy.optimize import brent, fmin
   #from scipy.special import erfc
   #if OccAvg is None:
      #nOcc = nElec/2
      #assert(nElec % 2 == 0)
      #return ORB_OCC * dot(Orbs[:,:nOcc],Orbs[:,:nOcc].T)
   #Type,beta = OccAvg
   #if Type == "Fermi":
      #MkOcc = lambda Mu: ORB_OCC/(1.+exp(beta*(Ew-Mu)))
      ## ^- fermi/dirac. Looks kinda complicated and has long tails.
   #elif Type == "Erfc":
      #MkOcc = lambda Mu: (.5*ORB_OCC)*erfc(beta*(Ew-Mu))
      ## ^- well.. it's symmetric. computationally nicer because it results in
      ##    less partially occupied orbitals.
   #elif Type == "Fixed":
      #Occ = zeros(Orbs.shape[0])
      #assert(sum(beta)==nElec)
      #Occ[:len(beta)] = beta
      #return mdot(Orbs,diag(Occ),Orbs.T)
   #else:
      #raise Exception("Occupation type not recognized. I know 'None', 'Fermi', and 'Erfc' ")
   #def Err(Mu):
      #return (sum(MkOcc(Mu)) - nElec)**2
   #Mu = brent(Err,brack=(Ew[nElec/2],Ew[nElec/2+1]),tol=1e-10)
   ## ^- note: The threshold needs to be real low! Otherwise DIIS blows up.
   #Occ = MkOcc(Mu)
   #Rdm = mdot(Orbs,diag(Occ),Orbs.T)
   #return Rdm


class FLatticeSystem(object):
   def __init__(self, WfDecl, SuperCell, InitialGuess, Params, Log):
      self.Model = SuperCell.Model
      self.SuperCell = SuperCell
      self.UnitCell = SuperCell.UnitCell
      assert(isinstance(Params, FMeanFieldParams))
      assert(isinstance(WfDecl, FWfDecl))
      self.P = Params
      self.WfDecl = WfDecl
      with Log.Section("lmf", "LATTICE HARTREE-FOCK", Log.SILENT):
         if self.WfDecl.OrbType == "UHF" and not self.SuperCell.SpinOrbs:
            raise Exception("unrestricted HF requires super-cells exposing a spin-orbital basis."
               "Please supply the OrbType='UHF' argument to FSuperCell.__init__ to make the super-cell"
               "compatible with UHF.")
         self.Init(InitialGuess, Log)

   def RunHf(self, Log):
      with Log.Section("lmf", "LATTICE HARTREE-FOCK",2):
         self._RunHf(Log)
      pass

   def Init(self, InitialGuess, Log):
      self.CoreEnergy = 0.
      # calculate the symmetry-unique elements of the hopping matrix,
      # in super-cell symmetrized real-space.
      sc = self.SuperCell
      self.CoreH = sc.MakeTijMatrix(sc.Sites, self.UnitCell)
      self.CoreH = self.CoreH.reshape(sc.nUnitCells, sc.nSitesU, sc.nSitesU)
      self.CoreH_Unpatched = 1. * self.CoreH

      self.Int2e_U = sc.MakeUiMatrix(sc.UnitCell)
      assert(self.Int2e_U.shape == (sc.nSitesU,))

      self.GuessType = 'CoreH'
      if InitialGuess is None:
         self.FockT = self.CoreH
      else:
         if isinstance(InitialGuess, FLatticeSystem):
            InitialGuess = InitialGuess.FockT
         if InitialGuess.shape == (sc.nSitesU,):
            # add a diagonal bias potential to the initial Fock matrix.
            # note that this can be used to break spatial symmetry, by
            # putting different potentials on alpha and beta electrons.
            self.FockT = 1. * self.CoreH
            self.FockT[0,:,:] += np.diag(InitialGuess)
            self.GuessType = 'CoreH+BiasPot'
         elif InitialGuess.shape == self.CoreH.shape:
            self.FockT = InitialGuess
            self.GuessType = 'Fock/Input'
         else:
            Log(Log.WARNING, "LatticeSystem::Init: Fock initial guess has wrong format. Using CoreH guess instead.")
            self.FockT = self.CoreH
      self.ORB_OCC = 1. if (self.WfDecl.OrbType == "UHF") else 2.

   def _EnergyFactor(self):
      #return 1./self.SuperCell.nClusterSites
      return self.SuperCell.EnergyFactor

   def _RunHf(self, Log):
      P = self.P
      sc = self.SuperCell
      Log(pInfoFmt,"Number of sites", sc.nSites)
      Log(pInfoFmt,"Number of sites in unit-cell", sc.nSitesU)
      Log()
      if self.P.MaxIt == 1:
         self.FockT = self.CoreH
         self.GuessType = "CoreH (due to non-SCF)"
      Log(pInfoFmt,"Initial guess", self.GuessType)
      f = 1./(1.*sc.nUnitCellsTotal)
      if self.WfDecl.OrbType == "RHF":
         Log(pInfoFmt,"Electrons in unit-cell", f*self.WfDecl.nElec)
      else:
         Log(pInfoFmt,"Electrons in unit-cell", (f*self.WfDecl.nElecA(), f*self.WfDecl.nElecB()) )
         #Log(pInfoFmt,"Electrons in unit-cell", self.WfDecl.nElecA()/(1.*sc.nUnitCellsTotal), "(alpha)")
         #Log(pInfoFmt,"Electrons in unit-cell", self.WfDecl.nElecB()/(1.*sc.nUnitCellsTotal), "(beta)")
      Log()

      # at this moment we have self.CoreH and self.Fock (initial guess).
      FockT = self.FockT
      Log("   {:^5s} {:^14s} {:^14s} {:^11s} {:>4s}",
          "ITER.","ENERGY","ENERGY CHANGE", "GRAD", "DIIS")
      EnergyFactor = self._EnergyFactor()
      # ^- calc energy per lattice model unit-cell
      LastRdm = None
      LastEnergy = 0.
      Converged = False
      dc = FDiisContext(P.DiisDim)
      for it in range(P.MaxIt):
         # transform Fock operator to translational irreps.
         FockK = sc.TransformToK(FockT)

         # make new orbitals and new rdm.
         Ew,Orbs = self.UpdateOrbitals(FockK)

         #Log("Fock spectrum (nSitesU x nGridK):\n{}", Ew.T)
         OccNumbers, RdmK, Mu, Gap = self.MakeRdm(Orbs, Ew)
         RdmT = sc.TransformToT(RdmK)

         # make new Fock matrix
         if P.MaxIt != 1 or True:
            FockT, j, k = self.MakeFock(RdmT, self.CoreH, self.Int2e_U, Orb=Orbs)
            self.HaveFock = True
         else:
            # in this case diagonalize CoreH only, do not incorporate
            # 2e potentials via HF mean-field, only via vcor.
            self.HaveFock = False
         Energy = (0.5 * (dot2(FockT, RdmT) + dot2(self.CoreH, RdmT)) + self.CoreEnergy).real
         Energy *= EnergyFactor

         #FockT = 1.*self.CoreH

         #print "Rdm:\n%s" % sc._ExpandFull(RdmT).real
         #print "Fock:\n%s" % sc._ExpandFull(FockT).real

         #Grd1 = np.dot(sc._ExpandFull(FockT),sc._ExpandFull(RdmT))
         #Grd1 = Grd1 - np.conj(Grd1.T)
         #print "Grd1:\n%s" % Grd1.real
         #raise SystemExit


         # calculate orbital gradient and energy
         FockK = sc.TransformToK(FockT)
         OrbGrad = np.zeros_like(FockK)
         for ik in xrange(FockK.shape[0]):
            OrbGrad[ik] = np.dot(FockK[ik],RdmK[ik])
            OrbGrad[ik] = OrbGrad[ik] - np.conj(OrbGrad[ik].T)
            # ^- fun fact: "OrbGrad -= OrbGrad.T" is an unsafe operation
            #    which, depending on circumstances, might or might not produce
            #    the expected result (well.. kinda obvious in hindsight).
            if self.WfDecl.OrbType == "UHF":
               OrbGrad[ik, ::2,1::2] = 0.
               OrbGrad[ik,1::2, ::2] = 0.

         fOrbGrad = la.norm(OrbGrad.flatten())
         dEnergy = Energy - LastEnergy
         LastEnergy = Energy

         Log(" {:5d} {:14.8f} {:+14.8f} {:11.2e} {:>6s}",
            (it+1 if P.MaxIt != 1 else 0), Energy, dEnergy, fOrbGrad, dc)

         if (fOrbGrad < P.ThrOrb and abs(dEnergy) < P.ThrDen):
            Converged = True
            break
         if (it == P.MaxIt - 1):
            break # not converged.

         # extrapolate
         #SkipDiis = (fOrbGrad > 1e-2 and it < 1)
         SkipDiis = it < P.DiisStart or fOrbGrad > P.DiisThr
         FockT, OrbGrad, c0 = dc.Apply(FockT, OrbGrad, Skip=SkipDiis)

      if (not Converged and P.MaxIt != 1):
         ErrMsg = "%s failed to converge."\
                  "  NIT =%4i  GRAD=%8.2e  DEN=%8.2e" % (self.WfDecl.OrbType, it+1, fOrbGrad, dEnergy)
         print "WARNING: %s" % ErrMsg
         if 0:
            print "Failed FSystem object stored to file 'brrrk' for investigation."
            import pickle
            pickle.dump(self,open('brrrk','w'))
         raise EHfConvergenceFail(ErrMsg)
      else:
         #print
         #PrintMatrix("Diagonal RDM:",diag(Rdm))
         Log()
         Log(pResultFmt, "Chemical potential", Mu)
         #Log(pResultFmt, "HOMO-LUMO gap", Gap)
         Log(pResultFmt, "Band gap", Gap)
         if Log['HfDensities']:
            Log()
            if self.WfDecl.OrbType == "UHF":
               Log(self.FmtRho,'Charge density', RdmT[0].real, '_')
               Log(self.FmtRho,'Spin density', RdmT[0].real, 'S')
            else:
               Log(self.FmtRho,'Charge density', RdmT[0].real, '_')
         Log()
         Log(pResultFmt, "1e energy", dot2(RdmT, self.CoreH)*EnergyFactor)
         Log(pResultFmt, "2e energy", .5*dot2(RdmT, FockT - self.CoreH)*EnergyFactor)
         Log(pResultFmt, "Hartree-Fock energy", Energy)

      Log(pResultFmt, "Total potential", FockT[0] - self.CoreH_Unpatched[0], self.WfDecl.OrbType)
      if P.MaxIt == 1:
         FockT = self.CoreH
         self.HaveFock = False
         # ^- other code assumes that .Fock is what we actually diagonalized.
         #Log(Log.WARNING, "Turning off Fock potentials in non-iterative lattice mean-field.")

      self.FockT = FockT
      self.RdmT = RdmT
      self.MeanFieldEnergy = EnergyFactor*Energy
      self.Mu = Mu
      self.Gap = Gap
      self.GuessType = 'Fock/Last'

   def FmtRho(self, Caption, Density, DensityType):
      return FmtRho(Caption, Density, DensityType, self.WfDecl.OrbType)
   def GetFock(self): return self.FockT
   def GetRdm(self): return self.RdmT

   def UpdateOrbitals(self, FockK):
      # format of output orbitals: nGridK x nSites x nOrb
      # i.e., MO dimension is last.
      sc = self.SuperCell
      assert(FockK.shape == (sc.nUnitCells, sc.nSitesU, sc.nSitesU))
      def UpdateOrbitalsInK(FockK):
         # diagonalize irreps in k-space individually
         Orbs = np.zeros_like(FockK)
         Ew = np.zeros((FockK.shape[0], Orbs.shape[2]))
         for ik in xrange(FockK.shape[0]):
            ew,ev = la.eigh(FockK[ik,:,:])
            Orbs[ik,:,:] = ev
            Ew[ik,:] = ew
         return Ew,Orbs
      #return UpdateOrbitalsInK(FockK)
      if self.WfDecl.OrbType == "RHF":
         return UpdateOrbitalsInK(FockK)
      else:
         ESC = ExtractSpinCompTK
         EwA,OrbsA = UpdateOrbitalsInK(ESC(FockK,0))
         EwB,OrbsB = UpdateOrbitalsInK(ESC(FockK,1))
         Ew = CombineSpinCompsTK(EwA,EwB)
         Orbs = CombineSpinCompsTK(OrbsA, OrbsB)
         return Ew,Orbs
      pass

   def MakeRdm(self, Orbs, Ews, nElec=None):
      sc = self.SuperCell

      if nElec is None:
         assert(Ews.shape == (sc.nUnitCells, sc.nSitesU))
         if self.WfDecl.OrbType == "UHF":
            # form RDM separately for alpha and beta orbitals.
            ESC = ExtractSpinCompTK
            OccNumbersA, RdmA, MuA, GapA = self.MakeRdm(ESC(Orbs,0), ESC(Ews,0), self.WfDecl.nElecA())
            OccNumbersB, RdmB, MuB, GapB = self.MakeRdm(ESC(Orbs,1), ESC(Ews,1), self.WfDecl.nElecB())
            CSC = CombineSpinCompsTK
            return CSC(OccNumbersA,OccNumbersB), CSC(RdmA,RdmB), (MuA,MuB), (GapA,GapB)
         else:
            nElec = self.WfDecl.nElec

      # find out which orbitals to occupy, and with which occupation numbers.
      OccNumbers = np.zeros_like(Ews)
      Ews_ = Ews.flatten()
      iEws = Ews_.argsort()

      # the number of orbitals we--in principle--would want to occupy.
      assert(nElec % self.ORB_OCC == 0)
      nOccFull = int(nElec / self.ORB_OCC + .5)

      Lumo = Ews_[iEws[nOccFull]]
      Homo = Ews_[iEws[nOccFull-1]] if (nOccFull != 0) else Lumo
      Mu = .5 * (Homo + Lumo)
      Gap = Lumo - Homo

      # first fill all orbitals below Mu-ThrDeg fully
      iOrbClo = Ews < Mu - self.P.ThrDeg
      iOrbOcc = Ews < Mu + self.P.ThrDeg
      nOrbClo = sum(iOrbClo.flatten())
      nOrbOcc = sum(iOrbOcc.flatten())
      nOrbAct = nOrbOcc - nOrbClo
      # now distribute the remaining electrons evenly amongst the degenerate
      # subspace
      nElecFull = self.ORB_OCC * nOrbClo
      nElecLeft = nElec - nElecFull
      assert(nElecLeft >= 0)
      if nOrbAct != 0:
         OccNumbers[iOrbOcc] = nElecLeft / nOrbAct
      OccNumbers[iOrbClo] = self.ORB_OCC

      if 0:
         print "Orbital energies (nSitesU x nGridK):\n{}".format(Ews.T - Mu)
         print "Occupation numbers (nSitesU x nGridK):\n{}".format(OccNumbers.T)
      #raise SystemExit

      #print "input orbitals:\n%s" % Orbs.T

      assert(np.all(OccNumbers >= 0.))
      # make the rdm: rdm[r,s] = \sum_i orb[r,i] occ-num[i] orb[s,i]
      #                        = \sum_i (orb[r,i] occ-num[i])**.5 * (orb[s,i] occ-num[i])**.5
      OrbsOcc = (Orbs.transpose((1,0,2)) * OccNumbers**.5).transpose((1,0,2))
      Rdm = np.zeros(Orbs.shape, OrbsOcc.dtype)
      for ik in xrange(Rdm.shape[0]):
         C = OrbsOcc[ik,:,:]   # <- nSites x nOrb
         Rdm[ik,:,:] = np.dot(C, np.conj(C.T))
         #print "\nC[ik={0}]:\n{1}\nRdm[ik={0}]:\n{2}".format(ik,C,Rdm[ik])
      #print "Rdm (nSitesU x nSitesU x nGridK):\n{}".format(Rdm.T)
      return OccNumbers, Rdm, Mu, Gap

   def MakeFock(self, RdmT, CoreH, Int2e_U, Orb=None):
      """make closed-shell/spin ensemble Fock matrix:
         f[rs] = h[rs] + j[rs] - .5 k[rs].
      where j[rs] = (rs|tu) rdm[tu]
      and   k[rs] = (rs|tu) rdm[su]
      Return (f,j,k).

      Args:
      #- RdmK: nGridK x nSitesU x nSitesU density matrix, in K-space
      - RdmT: nUnitCells x nSitesU x nSitesU density matrix, in T-space
      - CoreH: nUnitCells x nSitesU x nSitesU core Hamiltonian operator
        (i.e., sym-uq hopping matrix in real-space)
      - Int2e_U: nUnitCells x nSitesU vector of local 2e Coulomb on
        sites r.
      Notes:
      - This version supports local Hubbard-type coulomb interactions
        only. That means, in particular that
            (rs|tu) = U[r] * delta[r,s] * delta[r,t] * delta[r,u]
        and consequently k = j and both of them depend on the electron
        density rho[r] = rdm[r,r] only.
      - The Fock matrix is returned as (nUnitCells x nSitesU x nSitesU),
        i.e., with the first dimension unit-cell lattice translations
        in real-space.
      """
      # calculate density on sites in real-space.
      # Note that we only need the t=0 component for that,
      # as the interaction is local.
      sc = self.SuperCell
      if self.WfDecl.OrbType == "RHF":
         Rho = np.real(np.diag(RdmT[0,:,:])) # <- always real
         jk1 = np.diag(Int2e_U * Rho)
         jk = np.zeros_like(CoreH, Int2e_U.dtype)
         jk[0,:,:] = jk1
         return (CoreH + .5*jk, jk, jk)
      else:
         # calculate total and same-spin density on the sites.
         RhoSame = np.real(np.diag(RdmT[0,:,:]))
         RhoTotal = RhoSame[ ::2] + RhoSame[1::2]
         RhoTotal = np.hstack((RhoTotal,RhoTotal))
         #print "RhoSame:\n%s" % RhoSame
         #print "RhoTotal:\n%s" % RhoTotal
         #raise SystemExit
         j1 = np.diag(Int2e_U * RhoTotal)
         k1 = np.diag(Int2e_U * RhoSame)
         j = np.zeros_like(CoreH, Int2e_U.dtype); j[0,:,:] = j1
         k = np.zeros_like(CoreH, Int2e_U.dtype); k[0,:,:] = k1
         return (CoreH + (j-k), j, k)


FLatticeMeanField = FLatticeSystem

def flat01(Op):
   """transform (T,U,U) -> (T*U,U) matrix. (i.e., flatten dimensions 01)"""
   return Op.reshape( (Op.shape[0]*Op.shape[1],) + Op.shape[2:] )

def FmtRho(Caption, Density, DensityType, OrbType):
   (fa,fb) = {"A":(1.,0.), "B":(0.,1.), "_":(1.,1.), "S":(.5,-.5) }[DensityType]
   if len(Density.shape) == 2:
      Density = np.diag(Density)
   if OrbType == "UHF":
      A = Density[ ::2]
      B = Density[1::2]
   else:
      A = Density
      B = A
   R = fa*A + fb*B
   #print Density
   #raise SystemExit
   return ("{:36s}" + len(R)*" {:6.3f}").format(Caption, *R)




#def GetDefaultMeanfieldParameters():
   #Params = { "_lmf.ThrOrb": 1e-8,
              #"_lmf.ThrDen": 1e-8,
              ##"_lmf.MaxIt": 1024,
              #"_lmf.MaxIt": 12,
              #"_lmf.ThrDeg": 1e-5,
   #}
   #return Params


def _TestLatticeMeanFields():
   from numpy import set_printoptions, nan
   set_printoptions(precision=5,linewidth=10060,suppress=True,threshold=nan)

   from output import FOutputLog
   Log = FOutputLog()
   Params = { "t": 1.,
              "U": 1.,
              "_WfType": "RHF",
              "_lmf.WfDecl": FWfDecl(nElec=12),
              "_lmf.ThrOrb": 1e-8,
              "_lmf.ThrDen": 1e-8,
              #"_lmf.MaxIt": 1024,
              "_lmf.MaxIt": 12,
              "_lmf.ThrDeg": 1e-5,
   }
   Charge = 0
   if 1:
      L = 12
      nImp = 2
      Hub1d = FHubbardModel1d(1, Params)
      PhaseShift = [-1]
      #PhaseShift = [np.exp(.221j)]
      sc = FSuperCell(Hub1d, TotalSize=[L], PhaseShift=[-1], ClusterSize=[nImp])
   else:
      #FCls = FTestModel2d; Charge = 2

      FCls = FHubbardModel2dSquare; Charge = 0


      Lx,Ly=4,3
      #PhaseShift=[-1.,-1.]   # anti-periodic boundary conditions
      PhaseShift=[1.,1.]     # periodic boundary conditions
      #PhaseShift=[np.exp(.2j),np.exp(-.72j)] # wtf-like boundary conditions
      nImpX,nImpY = 2,1
      #nImpX,nImpY = 1,1
      nImp = nImpX * nImpY
      Hub2dSq = FCls((nImpX,nImpY), {'t': 1., 'U': 4.})
      sc = FSuperCell(Hub2dSq, [Lx/nImpX,Ly/nImpY], PhaseShift)

   Params["_lmf.WfDecl"] = FWfDecl(nElec=sc.nSites+Charge)
   lmf = FLatticeMeanField(sc, None, Params, Log)





if __name__ == "__main__":
   _TestLatticeMeanFields()



