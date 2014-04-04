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
import scipy.linalg as la
import itertools as it
from meanfield import FLatticeSystem, flat01
from quantum_systems import FSystem, FWfDecl
from helpers import mdot, MakeSmh, resmin, ExtractSpinComp, CombineSpinComps
from diis import FDiisContext
from output import *
from vcor_fit import FitVcorComponent

class FDmetParams:
   def __init__(self, THRBATHSVD=1e-6,
            THRDVC=1e-5, MAXIT=20, DIIS_START=2, DIIS_THR=1e-2, DIIS_DIM=4,
            VCOR_TYPE='Local', VCOR_FIT_TYPE='FullRdm',
            USE_INT2E=False):
      self.ThrBathSvd = THRBATHSVD # threshold in DMET SVD for inclusion of a bath orbital
      self.ThrDvc = THRDVC         # convergence threshold for change of correlation potential
      self.DiisStart = DIIS_START
      self.DiisDim = DIIS_DIM
      self.DiisThr = DIIS_THR
      self.VcorType = VCOR_TYPE
      self.VcorFitType = VCOR_FIT_TYPE
      self.MaxIt = MAXIT
      self.UseInt2e = USE_INT2E

class FDmetResult(object):
   def __init__(self, FragmentResults, ClusterFactor):
      # sum up the total energy and total electron number
      TotalEnergy = 0.
      TotalElec = 0.
      for (Fragment, r) in FragmentResults:
         Factor = Fragment.GetTotalFactor() * ClusterFactor
         TotalEnergy += Factor * r.EmbEnergy
         TotalElec += Factor * r.EmbElec
      self.TotalElecEnergy = TotalEnergy
      self.TotalEnergy = None
      self.TotalElec = TotalElec
      self.FragmentResults = FragmentResults
      self.FullSystemFock = None
      self.FullSystemVcor = None


class FReplica(object):
   """defines a replica of another fragment. If Trafo is not None,
   the other site's fragments are subjected to the transformation before
   being considered equivalent. Trafo must be unitary if given."""
   def __init__(self, Sites, Trafo=None):
      self.Sites = Sites
      if Trafo is None:
         self.Trafo = np.eye(len(self.Sites))
      else:
         self.Trafo = Trafo
         assert(la.norm((np.dot(Trafo, np.conj(Trafo.T)) - np.eye(len(self.Sites))).flatten()) < 1e-10)
   def __str__(self):
      return "[%s]" % ",".join(["%3i" % o for o in self.Sites])

class FFragment():
   def __init__(self, Sites, FullSystem, Method, Replicas=None):
      # Sites must be an integer sequence.
      self.Sites = np.array(Sites, int)
      self.Factor = 1.
      self.FullSystem = FullSystem
      self.Method = Method
      self.Replicas = [FReplica(self.Sites)]

   def MakeEmbeddedSystem(self, DmetParams):
      """construct the electronic system for the sites self.Sites as DMET
      embedded in the mean field result specified by FullSystemHfResult.
      Returns translation of self.Sites into the embedded system."""
      fs = self.FullSystem

      # make embedding basis and construct Hamiltonian
      # within it.
      EmbBasis = MakeEmbeddingBasis(self.Sites, flat01(fs.RdmT),
         fs.WfDecl.OrbType, DmetParams.ThrBathSvd)
      self.EmbBasis = EmbBasis

      # FIXME:
      #  - ToEmb() is actually really slow at this moment. And since it
      #    breaks translational invariance, this is to some degree intrinsic.
      #  - However, CoreH, CoreH_Unpatched, and Fock are all short-ranged
      #    and have support only on a few T vectors (and those are actually
      #    listed in SuperCell.MinSuperCell_iTs). For those, the transformation
      #    can be sped up by just ignoring all the other Ts.
      #  - The RDM does not have this property. However, we only transform
      #    it to find the number of electrons, and this can be done (although
      #    somewhat less savely) also without it.
      #    Additionally, we can just diagonalize the transformed F and will
      #    get the input RDM, as long as no bath sites are deleted and no
      #    degeneracies occur.

      # make CoreH, Int2e, and WfDecl
      CoreEnergy = 0   # (<- just added from full system)
      CoreH = self.ToEmb(fs.CoreH)

      nOrb = EmbBasis.shape[1]
      if not DmetParams.UseInt2e:
         # put local Us on the impurity sites and nothing on the rest
         Int2e_Ui = np.zeros(nOrb)
         Int2e_4ix = np.zeros((nOrb,nOrb,nOrb,nOrb))
         for (i,iSite) in enumerate(self.Sites):
            Int2e_Ui[i] = fs.Model.GetUi(fs.UnitCell[iSite])
            Int2e_4ix[i,i,i,i] = Int2e_Ui[i]
         if fs.WfDecl.OrbType == "UHF":
            assert(np.allclose(Int2e_Ui[ ::2]-Int2e_Ui[1::2],0.))
            Int2e_4ix[ ::2, ::2,1::2,1::2] = Int2e_4ix[ ::2, ::2, ::2, ::2]
            Int2e_4ix[1::2,1::2, ::2, ::2] = Int2e_4ix[ ::2, ::2, ::2, ::2]
      else:
         # transform U from full system basis into embedded basis
         Int2e_Ui_SuperCell = fs.SuperCell.MakeUiMatrix(fs.SuperCell.Sites)
         #print "Int2e_Ui_SuperCell.shape", Int2e_Ui_SuperCell.shape
         Int2e_4ix = MakeInt2e(Int2e_Ui_SuperCell, EmbBasis, EmbBasis)
         Int2e_Ui = None

      # add core Fock potential: That's the Fock potential of the electrons
      # not represented in the embedded system.
      #EmbRdm = self.ToEmb(flat01(fs.RdmT))
      EmbFock = self.ToEmb(flat01(fs.FockT))
      if self.FullSystem.HaveFock:
         CoreFock = EmbFock - self.MakeFock(EmbRdm, CoreH, Int2e_Ui,
            Int2e_4ix, OrbType=fs.WfDecl.OrbType)[0]
      else:
         CoreFock = 0.*EmbFock

      nImp = len(self.Sites)
      CoreH_Unpatched = self.ToEmb(fs.CoreH_Unpatched)
      CoreH[:nImp,:nImp] = CoreH_Unpatched[:nImp,:nImp]

      #def GetNumElec(Rdm):
         #fElec = np.trace(Rdm).real
         #nElec = int(fElec+.5)
         #assert(abs(nElec - fElec) < 1e-8)
         #return nElec
      def GetNumElec(nImp, Basis):
         nEmb = Basis.shape[1]
         nBath = nEmb - nImp
         return nBath * int(fs.WfDecl.OrbOcc())
         # ^- each bath orbital corresponds to one entangled occupied orbital
         #    and carries two electrons in the RHF case and one in the UHF case.
      if (fs.WfDecl.OrbType == "RHF" ):
         nElec = GetNumElec(nImp, EmbBasis)
         self.WfDecl = FWfDecl(nElec=nElec, OrbType="RHF")
      else:
         assert(fs.WfDecl.OrbType == "UHF")
         self.WfDecl = FWfDecl(nElecA=GetNumElec(nImp/2,ExtractSpinComp(EmbBasis,0)),
                               nElecB=GetNumElec(nImp/2,ExtractSpinComp(EmbBasis,1)),
                               OrbType=fs.WfDecl.OrbType)
      return FSystem(CoreEnergy, CoreH + CoreFock, Int2e_Frs=None, WfDecl=self.WfDecl,
                     Int2e_4ix=Int2e_4ix,
                     CoreFockV=CoreFock, FockGuess=EmbFock),\
             list(range(len(self.Sites))) # <- first sites in embedded sys are the impurity sites

   def MakeFock(self, Rdm, CoreH, Int2e_Ui, Int2e_4ix=None, Orb=None, OrbType=None):
      if OrbType is None:
         OrbType = self.WfDecl.OrbType
      if Int2e_Ui is not None:
         # see meanfield.py.
         if OrbType == "RHF":
            Rho = np.real(np.diag(Rdm)) # <- always real
            jk = np.diag(Int2e_Ui * Rho)
            return (CoreH + .5*jk, jk, jk)
         else:
            # calculate total and same-spin density on the sites.
            RhoSame = np.real(np.diag(Rdm))
            RhoTotal = RhoSame[ ::2] + RhoSame[1::2]
            RhoTotal = np.hstack((RhoTotal,RhoTotal))
            j = np.diag(Int2e_Ui * RhoTotal)
            k = np.diag(Int2e_Ui * RhoSame)
            return (CoreH + (j-k), j, k)
      elif Int2e_4ix is not None and OrbType == "RHF":
         j = np.einsum('rstu,tu', Int2e_4ix, Rdm)
         k = np.einsum('rtsu,tu', Int2e_4ix, Rdm)
         return (CoreH + j -.5*k, j, k)
      else:
         assert(0)

   def ToEmb(self, M):
      """Transform matrix from full system into the embedded system."""
      #return mdot(self.EmbBasis.T, M, self.EmbBasis)
      return self.FullSystem.SuperCell.TransformToEmb(self.EmbBasis, M, self.EmbBasis)

   def Run(self, DmetParams, Log):
      FragmentSystem, EmbSites = self.MakeEmbeddedSystem(DmetParams)
      Result = FragmentSystem.Run(self.Method, SpecialSites=EmbSites, Log=Log)
      return Result

   def FitCorrelationPotential(self, EmbSystemHlResult, FullSystemHfResult,
         DmetParams, OrbType, Log):
      """fit a matrix nEmb x nEmb, referring to an operator connecting
      all embedded sites to each other, to more closely approach the
      correlated calculation result in with the mean-field method.

      Returns nEmb x nEmb matrix (always), even if due to the fitting
      criterion used the actual matrix is only a subset of this."""
      # transform Fock matrix from full system into embedding basis.
      g = DmetParams
      EmbFock = self.ToEmb(FullSystemHfResult.GetFock())
      nImp = len(self.Sites)
      nEmb = self.EmbBasis.shape[1]
      RdmHl = EmbSystemHlResult.GetRdm()  # <- high level RDM we intent to fit.
      assert(self.WfDecl.OrbType == OrbType)
      ORB_OCC = self.WfDecl.OrbOcc()

      if ( self.WfDecl.OrbType == "RHF" ):
         assert(self.WfDecl.nElec%2 == 0)
         return FitVcorComponent(EmbFock, nImp, (1./ORB_OCC)*RdmHl, g.VcorType, g.VcorFitType)
      else:
         assert(self.WfDecl.OrbType == "UHF")
         return CombineSpinComps(\
            FitVcorComponent(EmbFock[ ::2, ::2], nImp/2, RdmHl[ ::2, ::2], g.VcorType, g.VcorFitType),
            FitVcorComponent(EmbFock[1::2,1::2], nImp/2, RdmHl[1::2,1::2], g.VcorType, g.VcorFitType))
      pass

   def GetTotalFactor(self):
      return self.Factor * len(self.Replicas)

def MakeInt2e(Ui, EmbBasisIJ, EmbBasisKL):
    """transform local on-site interaction from large lattice into embedding basis.
    Result is a 4-tensor (ij|kl) containing the matrix elements of the 2e Hamiltonian
        (ij|kl) E^i_j E^k_l
    in the embedding basis.
    EmbBasisIJ contains the basis used to transform the (ij| pairs, EmbBasisKL contains
    # the basis for the |kl) pairs.
    """
    assert(EmbBasisIJ.shape == EmbBasisKL.shape)
    assert(Ui.shape == (EmbBasisIJ.shape[0],))
    nOrb = EmbBasisIJ.shape[1]
    Int2e = np.zeros( (nOrb,nOrb,nOrb,nOrb) )
    def outer4(a,b,c,d):
       na = np.newaxis
       return ((a[:,na] * b[na,:])[:,:,na] * c[na,:])[:,:,:,na] * d[na,:]
    for iSite in range(EmbBasisIJ.shape[0]):
       u = EmbBasisIJ[iSite,:]
       v = EmbBasisKL[iSite,:]
       # see nsite.py. This sums U[ijkl] := \sum_mu U[mu]*c[mu,i]*c[mu,j]*c[mu,k]*c[mu,l].
       Int2e[:,:,:,:] += Ui[iSite]*outer4(u,u,v,v)
    return Int2e

def MakeFockX1(Rdm, CoreH, Int2e_Ui, Int2e_4ix=None, Orb=None, OrbType="RHF"):
   if Int2e_Ui is not None:
      # see meanfield.py.
      if OrbType == "RHF":
         Rho = np.real(np.diag(Rdm)) # <- always real
         jk = np.diag(Int2e_Ui * Rho)
         return (CoreH + .5*jk, jk, jk)
      else:
         # calculate total and same-spin density on the sites.
         RhoSame = np.real(np.diag(RdmT[0,:,:]))
         RhoTotal = RhoSame[ ::2] + RhoSame[1::2]
         RhoTotal = np.hstack((RhoTotal,RhoTotal))
         j1 = np.diag(Int2e_U * RhoTotal)
         k1 = np.diag(Int2e_U * RhoSame)
         return (CoreH + (j-k), j, k)
   elif Int2e_4ix is not None and OrbType == "RHF":
      j = np.einsum('rstu,tu', Int2e_4ix, Rdm)
      k = np.einsum('rtsu,tu', Int2e_4ix, Rdm)
      return (CoreH + j -.5*k, j, k)
   else:
      assert(0)

def MakeEmbeddingBasis1(Sites, Rdm, ThrBathSvd):
   """calculate the DMET embedding basis for the given list of sites (integer
   array). Rdm must be an idempotent density matrix (that is, a density matrix
   from a 1e matrix diagonalization)."""
   if (type(Sites) is not np.ndarray):
      Sites = array(Sites, int)

   # make proto-bath vectors.
   vecs0 = 1.*Rdm[:,Sites]
   vecs0[Sites,:] = 0
   # orthogonalize
   S = np.dot(vecs0.T, vecs0)
   if 0:
      vecs0 = np.dot(vecs0, MakeSmh(S))
   else:
      ew,ev = la.eigh(-S)
      nBath = len([o for o in ew if -o > ThrBathSvd])
      if nBath == len(Sites):
         vecs0 = np.dot(vecs0, MakeSmh(S))
      else:
         #raise Exception("Had to delete bath functions.... Smallest eigenvalue: %8.2e" % ew[-1])
         vecs0 = np.dot(vecs0, np.dot(ev[:,:nBath],np.diag((-ew[:nBath])**-.5)))

   nImp = len(Sites)
   basis = np.zeros((Rdm.shape[0], nImp + vecs0.shape[1]), Rdm.dtype)
   basis[Sites,:nImp] = np.eye(nImp)
   basis[:,nImp:] = vecs0
   return basis


def MakeEmbeddingBasis(Sites, Rdm, OrbType, ThrBathSvd):
   if OrbType == "RHF":
      return MakeEmbeddingBasis1(Sites, Rdm, ThrBathSvd)
   else:
      # FIXME: hmpf... might have different number of bath sites for alpha
      #        and beta... wat do?
      basisA = MakeEmbeddingBasis1(np.array(Sites[ ::2])/2, ExtractSpinComp(Rdm,0), ThrBathSvd)
      basisB = MakeEmbeddingBasis1(np.array(Sites[1::2])/2, ExtractSpinComp(Rdm,1), ThrBathSvd)
      assert(basisA.shape == basisB.shape)
      return CombineSpinComps(basisA, basisB)

#class F

#class FDmetResult:
   


#class FDmetContext(object):
   #def __init__(self, LatticeSystem, Fragments, DmetParams):
      #self.LatticeSystem = LatticeSystem
      #self.Fragments = []
      #for F in Fragments:
         #self.Fragments.append( FFragment(Method=F[0],
            #Sites=F[1], FullSystem=self.LatticeSystem, Replicas=None) )
      #self.DmetParams = DmetParams

   #def Run(self):
      #LatticeSystem.RunHf(Log=Log)

      ## fixme: should not be set here.
      #FragmentSites = list(range(self.Params.SuperCell.ClusterSize[0]))
      #Fragment = FFragment(Sites=FragmentSites, FullSystem=LatticeSystem, Method='FCI')

      #FragmentSystem, EmbSites = Fragment.MakeEmbeddedSystem(self.Params._DMET)
      #FragmentSystem.RunHf(Log, SpecialSites=Fragment.Sites)
      #FragmentSystem.RunFci(Log, SpecialSites=Fragment.Sites)


class FDmetContext:
   def __init__(self, LatticeSystem, Fragments, DmetParams, Log):
      """represents the entire system and its subdivision into fragments."""
      assert(isinstance(LatticeSystem, FLatticeSystem))
      assert(isinstance(DmetParams, FDmetParams))
      self.FullSystem = LatticeSystem
      self.DmetParams = DmetParams
      self.Fragments = []
      def AdaptSiteList(L):
         if ( LatticeSystem.SuperCell.SpinOrbs and not
              LatticeSystem.Model.SitesAreSpinOrbitals ):
            # input model is a spatial model, but we're doing UHF. Duplicate
            # all input sites to get both the alpha orbitals at the even indices
            # and the beta orbitals at the odd indices.
            return list(it.chain(*[(2*o, 2*o+1) for o in L]))
         else:
            # we either have spatial sites, or the model is explicitly a spin-
            # orbital model, in which case probably the input fragmentation is
            # also given in terms of spin-orbitals already. In either case no
            # adjustment is required.
            return L
      for F in Fragments:
         self.Fragments.append( FFragment(Method=F[0],
            Sites=AdaptSiteList(F[1]), FullSystem=self.FullSystem,
            Replicas=None) )
      self.CheckFragmentation(Log)

   def CheckFragmentation(self, Log):
      # check if each site is covered exactly once.
      SiteOccurence = {}
      MaxSite = -1
      for Fragment in self.Fragments:
         for Replica in Fragment.Replicas:
            for Site in Replica.Sites:
               if Site not in SiteOccurence:
                  SiteOccurence[Site] = 0.
               SiteOccurence[Site] += Fragment.Factor
               if Site > MaxSite: MaxSite = Site
      Ok = True
      if ( MaxSite == -1 ):
         Log(Log.WARNING, "Fragmenation covers no sites.")
         Ok = False
      if ( MaxSite + 1 != len(self.FullSystem.UnitCell) ):
         Log(Log.WARNING, "Some sites of unit cell not present in fragmentation:"+\
               "\n         nSitesU = {}  MaxSite = {}.", len(self.FullSystem.UnitCell), MaxSite)
         Ok = False
      for Site in range(MaxSite):
         if Site not in SiteOccurence:
            Log(Log.WARNING, "Site {} not covered in fragmentation!", Site)
            Ok = False
         if abs(SiteOccurence[Site] - 1.) > 1e-8:
            Log(Log.WARNING, "Site {} is covered {:.5f} times. Factor sum should be unity!", Site, CoveredSites[Site])
            Ok = False
      if not Ok:
         Msg = "There are problems with the fragmentation. Check and ignore."
         raise Exception(Msg)
         Log(Log.ERROR, Msg)

   def RunFragments(self, FullSystemHfResult, Log):
      """run HL electronic structure calculation on the individual
      fragments. return results of each of them."""
      FragmentResults = []
      PrintEmb = Log['FragmentCalc']
      for (iFragment,Fragment) in enumerate(self.Fragments):
         #if PrintEmb:
            #Log(FmtHeading("EMBEDDED SUBSYSTEM: Id =%3i  F =%9.5f"\
               #% (iFragment, Fragment.Factor), 1))
         with Log.Section('F%02x'%iFragment, "EMBEDDED SUBSYSTEM: Id =%3i  F =%9.5f" % (iFragment, Fragment.Factor), 1):
            sr = Fragment.Run(self.DmetParams, Log)
            FragmentResults.append((Fragment, sr))

      ClusterFactor = self.FullSystem._EnergyFactor()
      # ^- for normalization with super-cell cluster size. Purely cosmetic.
      return FDmetResult(FragmentResults, ClusterFactor)

   def FitCorrelationPotentials(self, FragmentGroupResult, FullSystemHfResult, Log):
      FragmentPotentials = []
      for (Fragment,r) in FragmentGroupResult.FragmentResults:
         vloc = Fragment.FitCorrelationPotential(r, FullSystemHfResult,
            self.DmetParams, self.FullSystem.WfDecl.OrbType, Log=Log)
         FragmentPotentials.append(vloc)
      FragmentGroupResult.FragmentPotentials = FragmentPotentials

   #def PrintFragmentResults(self, FragmentGroupResult, iMacroIter, dVsum, PrintEmb):
      #"""print summary of the fragment results."""
      #jgr = FragmentGroupResult
      #assert(isinstance(jgr,FDmetResult))

      #PrintHeading("FRAGMENT SUMMARY & SYSTEM REASSEMBLY -- ITER. %i"\
         #% (1+iMacroIter), int(PrintEmb))
      #print ("  {:^5s} {:^13s} {:^14s}{:^12s}   {:>6s}").format(
         #"SUB.SYS","FACTOR","ELEC.ENERGY", "ELECTRONS", "d[V]")
      #for (iFragment, (Fragment,r)) in enumerate(jgr.FragmentResults):
         #dV = "  --"
         #if ( jgr.FragmentPotentials is not None ):
            #dV = "%9.2e" % norm(jgr.FragmentPotentials[iFragment])
         #print " %5i %13.5f %15.8f %11.6f %10s" % (iFragment,
            #Fragment.Factor, r.EmbEnergy, r.EmbElec, dV)
      #print
      #print g.pResultFmt % ("Potential convergence", dVsum)
      #print
      #print g.pResultFmt % ("Number of electrons", jgr.TotalElec)
      #print g.pResultFmt % ("Total electronic energy", jgr.TotalElecEnergy)
      #print g.pResultFmt % ("Full.sys. core energy", self.FullSystem.CoreEnergy)
      #print g.pResultFmt % ("Total energy", jgr.TotalEnergy)
      #print

   def Run(self, Log, StartingGuess={}):
      PrintEmb = Log['EmbeddedCalc']

      FullCoreH_Unpatched = self.FullSystem.CoreH_Unpatched
      if 'vcor' in StartingGuess:
         VcorLarge = 1.*StartingGuess['vcor']
      else:
         VcorLarge = np.zeros_like(FullCoreH_Unpatched[0,:,:])
      dc = FDiisContext(self.DmetParams.DiisDim)

      g = self.DmetParams

      def RunVcorUpdate(VcorLarge):
         #assert(np.allclose(VcorLarge[1::2, ::2],0.))
         #assert(np.allclose(VcorLarge[ ::2,1::2],0.))
         with Log.Section("lmf","FULL SYSTEM -- MEAN FIELD CALCULATION", 1):
            #vFock = self.FullSystem.FockT[0] - self.FullSystem.CoreH[0]
            self.FullSystem.CoreH[0] = FullCoreH_Unpatched[0] + VcorLarge
            #self.FullSystem.FockT[0] = self.FullSystem.CoreH[0] + vFock
            self.FullSystem.FockT[0] = FullCoreH_Unpatched[0] + VcorLarge
            # ^- this adds the changed vloc also to the FockT initial guess.
            #    required if iterations are disabled.
            FullSystemHfResult = self.FullSystem._RunHf(Log=Log)
            FullSystemHfResult = self.FullSystem

         jgr = self.RunFragments(FullSystemHfResult, Log)
         jgr.MeanFieldEnergy = FullSystemHfResult.MeanFieldEnergy
         jgr.TotalEnergy = jgr.TotalElecEnergy + self.FullSystem.CoreEnergy
         self.FitCorrelationPotentials(jgr, FullSystemHfResult, PrintEmb)

         # now: add changed potentials to large system, and run again.
         dVcorLarge = np.zeros(VcorLarge.shape)
         assert(g.VcorType == "Local" or g.VcorType == "Diagonal")
         for (Fragment, VcorFull) in zip(self.Fragments, jgr.FragmentPotentials):
            nImp = len(Fragment.Sites)
            Vcor = Fragment.Factor * VcorFull[:nImp,:nImp]
            for Replica in Fragment.Replicas:
               VcorR = mdot(Replica.Trafo, Vcor, Replica.Trafo.T)
               E = list(enumerate(Replica.Sites))
               for (i,iSite),(j,jSite) in it.product(E,E):
                  dVcorLarge[iSite,jSite] += VcorR[i,j]

         return dVcorLarge, jgr, FullSystemHfResult

      with Log.Section('vc', "SELF-CONSISTENT EMBEDDING: CORR.POTENTIAL FIT", 1):
         IterationHistory = []
         IterationHistory.append("   THRDVC = {0:8.2e}       VCOR_TYPE = {1}/{2}\n".format(g.ThrDvc, g.VcorType, g.VcorFitType))
         IterationHistory.append("   {:^5s} {:^15s} {:^14s}  {:>6s}  {:>6s}".format("ITER.","TOT.ENERGY", "ELECTRONS", "d[V]", "DIIS"))
         for iMacroIter in range(g.MaxIt):
            with Log.Section("%02x" % iMacroIter, "CORRELATION POTENTIAL FIT - ITER. %i" % iMacroIter):
               dVcorLarge, jgr, FullSystemHfResult = RunVcorUpdate(VcorLarge)
               dVsum = la.norm(dVcorLarge.flatten())
               if (iMacroIter == 0 ):
                  InitialHfEnergy = FullSystemHfResult.MeanFieldEnergy

               DiisText = " -  -"
               if dc: DiisText = "%2i %2i" % (dc.nDim, dc.iNext)
               IterationHistory.append(" %5i %15.8f %13.7f    %8.2e  %s" % (iMacroIter+1,
                  jgr.TotalEnergy, jgr.TotalElec, dVsum, DiisText))

               #self.PrintFragmentResults(jgr, iMacroIter, dVsum, PrintEmb)

               jgr.dVc = dVsum
               if ( dVsum < g.ThrDvc ):
                  Log("Convergence criteria met -- done.")
                  break

               Log("PROGESS OF CORR.POTENTIAL FIT:\n")
               for Line in IterationHistory:
                  Log(Line)


               vFock = np.zeros((1))
               vFock = self.FullSystem.FockT[0] - self.FullSystem.CoreH[0]
               SkipDiis = dVsum > self.DmetParams.DiisThr and iMacroIter < self.DmetParams.DiisStart
               print "diis information"
               print VcorLarge
               print dVcorLarge
               print vFock
               VcorLarge, dVcorLarge, vFock, c0 = dc.Apply(VcorLarge, dVcorLarge, vFock, Skip=SkipDiis)
               print "diis information"
               print VcorLarge
               print dVcorLarge
               print vFock
               if not SkipDiis:
                  Log("\nVcor extrapolation: DIIS{{{:2d} {:2d}  {:8.2e}}}", dc.nDim, dc.iNext, c0)
               self.FullSystem.FockT[0] = vFock + self.FullSystem.CoreH[0]
               #self.FullSystem.FockT = Fock
               VcorLarge += dVcorLarge
               #raise Exception("absicht.")
         Log("")
         Log("DMET SCF cycle converged after {} iterations", 1+iMacroIter)
         Log("Final residual on unit-cell vcor is {:.2e}", dVsum)
         #Log("Fit history:")
         Log()
         for Line in IterationHistory:
            Log(Line)
         Log()
         DiagOnly = False
         if DiagOnly:
            DiagDesc,DiagMap = " (diagonal)", lambda x:np.diag(x)
         else:
            DiagDesc,DiagMap = "", lambda x:x
         Log(pResultFmt, "Final correlation potential%s"%DiagDesc, DiagMap(VcorLarge), self.FullSystem.WfDecl.OrbType)
         Log(pResultFmt, "Final mean field rdm%s"%DiagDesc, DiagMap(FullSystemHfResult.RdmT[0].real), self.FullSystem.WfDecl.OrbType)
         Log()
         Log(pResultFmt, "Mean-field energy", InitialHfEnergy, "initial")
         Log(pResultFmt, "Mean-field energy", FullSystemHfResult.MeanFieldEnergy, "final")
         Log()
         Log(pResultFmt, "DMET energy", jgr.TotalEnergy)
         #self.PrintFragmentResults(jgr, iMacroIter, dVsum, PrintEmb)

      jgr.FullSystemFock = FullSystemHfResult.GetFock()
      jgr.FullSystemVcor = VcorLarge
      return jgr




