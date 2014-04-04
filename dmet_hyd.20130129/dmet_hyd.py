# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Gerald Knizia and Garnet K.-L. Chan, 2012;
# (c) Princeton University, 2012

import itertools as it
itt = it
from os import path
from sys import argv
from copy import copy, deepcopy
from numpy import *
from numpy.random import random as uniform_random, seed
from scipy.linalg import *

from diis import FDiisContext
from helpers import *
import settings as g
from quantum_systems import *
from molecules import *
from vcor_fit import FitVcorComponent

seed(6544)

STOP_AFTER_FRAGMENT = None

class FReplica(object):
   def __init__(self, Sites, Trafo=None):
      self.Sites = Sites
      if Trafo is None:
         self.Trafo = eye(len(self.Sites))
      else:
         self.Trafo = Trafo
         assert(norm((dot(Trafo, Trafo.T) - eye(len(self.Sites))).flatten()) < 1e-10)
   def __str__(self):
      return "[%s]" % ",".join(["%3i" % o for o in self.Sites])

class FFragmentDesc(object):
   def __init__(self, Sites, Factor, Replicas=None):
      """defines fragment in terms of sites of the entire system.

      If Replicas is provided, it must be a list of lists of sites,
      all of the same length, which are understood as being identical
      to the main system (defined by Sites). Correlation potentials
      will be replicated to those sites, and the total energy and
      electron number will be .Factor * len(Replicas). Must include
      .Sites itself."""
      def SplitSpinComponents(L):
         if g.WF_TYPE == "RHF":
            # do nothing.
            return L
         else:
            assert(g.WF_TYPE == "UHF")
            # add the spatial site's alpha (even i) and beta (odd i)
            # spin-orbital
            return list(it.chain(*[(2*o,2*o+1) for o in L]))

      self.Sites = SplitSpinComponents(Sites)
      self.Factor = Factor
      if Replicas is None:
         self.Replicas = [FReplica(self.Sites)]
      else:
         # if input is a site list (i.e., no trafo), make replica objects
         # out of it.
         def ToReplica(L):
            if isinstance(L,FReplica): return L
            else: return FReplica(L)
         self.Replicas = map(ToReplica, Replicas)
         self.Replicas = map(SplitSpinComponents, self.Replicas)
         nSites = len(self.Replicas[0].Sites)
         for Replica in self.Replicas:
            if ( len(Replica.Sites) != nSites ):
               raise Exception("Replicas don't all have the same length. nUniqueBaseFragments correct?")
         if g.WF_TYPE == "UHF":
            print "\nWARNING: inverting site order of every second replica in UHF mode!"
            for i in range(len(self.Replicas)):
               if i % 2 == 0 and g.WF_TYPE == "UHF":
                  self.Replicas.Sites[i] = self.Replicas.Sites[i][::-1]
            print "         new replicas: %s\n" % self.Replicas
         iSelf = [tuple(o.Sites) for o in self.Replicas].index(tuple(self.Sites))
         r0 = 1. * self.Replicas[iSelf].Trafo
         for r in self.Replicas:
            r.Trafo = dot(r.Trafo,r0.T)
         #assert(all(self.Replicas[iSelf].Trafo == eye(len(self.Sites))))

   def GetTotalFactor(self):
      return self.Factor * len(self.Replicas)

def PrintSiteCoverage(Basis, BasisDesc):
   dummy = zeros(Basis.shape[0])
   for i in range(Basis.shape[1]):
      dummy += Basis[:,i]**2
   dummy = sqrt(dummy)
   print "site coverage of the embedding basis:\n"
   print "  Basis:    %s" % " ".join(["  %4i  " % i for o in range(Basis.shape[1])])
   for i in range(len(dummy)):
      Individual = " ".join(["%8.4f" % o**2 for o in Basis[i,:]])
      print "%14s: %s  %8.4f" % (BasisDesc[i], Individual, dummy[i])
   print
   print "BasisLabels = [%s]" % ", ".join(['"%s"' % o for o in BasisDesc])
   print "Full embedding basis:\n" + str(Basis)
   print


def MakeEmbeddingBasis1(Sites, rdm):
   """calculate the DMET embedding basis for the given list of sites (integer
   array). rdm must be an idempotent density matrix (that is, a density matrix
   from a 1e matrix diagonalization)."""
   if (type(Sites) is not ndarray):
      Sites = array(Sites, int)
   # make proto-bath vectors.
   vecs0 = 1.*rdm[:,Sites]
   vecs0[Sites,:] = 0
   # orthogonalize
   S = dot(vecs0.T, vecs0)
   if 0:
      ew,ev = eigh(-S); ew *= -1
      print "       bath ew: [%s]" % " ".join("%12.8f" % o for o in ew)
   if 0:
      vecs0 = dot(vecs0, MakeSmh(S))
   else:
      ew,ev = eigh(-S)
      nBath = len([o for o in ew if -o > 1e-8])
      if nBath == len(Sites):
         vecs0 = dot(vecs0, MakeSmh(S))
      else:
         #raise Exception("Had to delete bath functions.... Smallest eigenvalue: %8.2e" % ew[-1])
         vecs0 = dot(vecs0, dot(ev[:,:nBath],diag((-ew[:nBath])**-.5)))

   nImp = len(Sites)
   basis = zeros((rdm.shape[0], nImp + vecs0.shape[1]))
   basis[Sites,:nImp] = eye(nImp)
   basis[:,nImp:] = vecs0

   if 0:
      dummy = zeros(rdm.shape[0])
      for i in range(basis.shape[1]):
         dummy += basis[:,i]**2
      dummy = sqrt(dummy)
      print "site coverage of the embedding basis:\n%s" % dummy
   return basis


def MakeEmbeddingBasis(Sites, rdm):
   if g.WF_TYPE == "RHF":
      return MakeEmbeddingBasis1(Sites, rdm)
   else:
      # FIXME: hmpf... might have different number of bath sites for alpha
      #        and beta... wat do?
      basisA = MakeEmbeddingBasis1(array(Sites[ ::2])/2, ExtractSpinComp(rdm,0))
      basisB = MakeEmbeddingBasis1(array(Sites[1::2])/2, ExtractSpinComp(rdm,1))
      return CombineSpinComps(basisA, basisB)


class FFragment(FFragmentDesc, FSystem):
   def __init__(self, Sites, Factor, FullSystem, Method, Replicas=None):
      super(FFragment,self).__init__(Sites, Factor, Replicas)
      self.FullSystem = FullSystem
      self.Method = Method

   def MakeEmbeddedSystem(self, FullSystemHfResult):
      """construct the electronic system for the sites self.Sites as DMET
      embedded in the mean field result specified by FullSystemHfResult.
      Returns translation of self.Sites into the embedded system."""
      fs = self.FullSystem
      fr = FullSystemHfResult

      # make embedding basis and construct Hamiltonian
      # within it.
      EmbBasis = MakeEmbeddingBasis(self.Sites, fr.GetRdm())
      self.EmbBasis = EmbBasis
      if 0:
         Fock = self.ToEmb(fr.GetFock())
         Rdm = self.ToEmb(fr.GetRdm())
         print mdot(Fock,Rdm) - mdot(Rdm,Fock)
         raise SystemExit

      # make CoreH, Int2e, and WfDecl
      CoreEnergy = 0   # (<- just added from full system)
      CoreH = self.ToEmb(fs.CoreH)
      #Int2e_Frs = einsum("Frs,ri,sj", fs.Int2e_Frs, EmbBasis, EmbBasis)
      # ^- doesn't seem like einsum bothers with factorizing this...
      Frj = einsum("Frs,sj->Frj", fs.Int2e_Frs, EmbBasis)
      Int2e_Frs = einsum("ri,Frj->Fij", EmbBasis, Frj)

      nImp = len(self.Sites)
      if not g.USE_INT2E:
         Int2e_Frs[:,nImp:,:] = 0
         Int2e_Frs[:,:,nImp:] = 0

      # add core Fock potential: That's the Fock potential of the electrons
      # not represented in the embedded system.
      EmbRdm = self.ToEmb(fr.GetRdm())
      EmbFock = self.ToEmb(fr.GetFock())
      CoreFock = EmbFock - MakeFockOp(EmbRdm, CoreH, Int2e_Frs)[0]
      if STOP_AFTER_FRAGMENT is not None:
         PrintSiteCoverage(EmbBasis, self.FullSystem.BasisDesc)
         print "Input EmbRdm:\n%s\nInput EmbFock:\n%s" % (EmbRdm, self.ToEmb(fr.GetFock()))
      #print "CoreH:\n%s\nCoreFock:\n%s\n" % (CoreH, CoreFock)

      if 1:
         p = self.Sites[0] == 0
         CoreH_Unpatched = self.ToEmb(fs.CoreH_Unpatched)
         #if p: print "DIFF -- COREH - COREH_Unpatched:\n%s" % (CoreH - CoreH_Unpatched)
         CoreH[:nImp,:nImp] = CoreH_Unpatched[:nImp,:nImp]

      def GetNumElec(Rdm):
         fElec = trace(Rdm)
         nElec = int(fElec+.5)
         assert(abs(nElec - fElec) < 1e-8)
         return nElec
      if (g.WF_TYPE == "RHF" ):
         nElec = GetNumElec(EmbRdm)
         #assert(nElec == int(g.ORB_OCC)*len(self.Sites))
         # ^- this holds as long as no bath orbitals are deleted.
         self.WfDecl = FWfDecl(nElec=nElec)
      else:
         assert(g.WF_TYPE == "UHF")
         self.WfDecl = FWfDecl(nElecA=GetNumElec(ExtractSpinComp(EmbRdm,0)),
                               nElecB=GetNumElec(ExtractSpinComp(EmbRdm,1)))
      return FSystem(CoreEnergy, CoreH + CoreFock, Int2e_Frs, self.WfDecl,
                     CoreFockV=CoreFock, FockGuess=EmbFock),\
             list(range(len(self.Sites))) # <- first sites in embedded sys are the impurity sites

   def ToEmb(self, M):
      """Transform matrix from full system into the embedded system."""
      return mdot(self.EmbBasis.T, M, self.EmbBasis)

   def Run(self, FullSystemHfResult, Print = True):
      FragmentSystem, EmbSites = self.MakeEmbeddedSystem(FullSystemHfResult)
      Result = FragmentSystem.Run(self.Method, SpecialSites = EmbSites, Print = Print)
      Result.EmbSystem = FragmentSystem # <- FIXME: just for experiments. Remove.
      return Result

   def FitCorrelationPotential(self, EmbSystemHlResult, FullSystemHfResult, Print = True):
      """fit a matrix nEmb x nEmb, referring to an operator connecting
      all embedded sites to each other, to more closely approach the
      correlated calculation result in with the mean-field method.

      Returns nEmb x nEmb matrix (always), even if due to the fitting
      criterion used the actual matrix is only a subset of this."""

      # transform Fock matrix from full system into embedding basis.
      EmbFock = self.ToEmb(FullSystemHfResult.GetFock())
      nImp = len(self.Sites)
      nEmb = self.EmbBasis.shape[1]
      RdmHl = EmbSystemHlResult.GetRdm()  # <- high level RDM we intent to fit.
      nOcc = self.WfDecl.nOcc

      def FitComponent(EmbFock, nImp, nEmb, RdmHl, nOcc):
         assert(g.VLOC_TYPE in ["Local","Diagonal"])
         vloc = zeros((nEmb,nEmb))
         vloc[:,:] = FitVcorComponent(EmbFock, nImp, RdmHl, g.VLOC_TYPE, g.VLOC_FIT_TYPE)
         return vloc

      if ( g.WF_TYPE == "RHF" ):
         return FitComponent(EmbFock, nImp, nEmb, RdmHl/2, nOcc)
         # ^- RdmHl/2: input is supposed to be approximately idempotent, i.e.,
         #    normalized to ORB_OCC==1..
      else:
         assert(g.WF_TYPE == "UHF")
         return CombineSpinComps(\
            FitComponent(EmbFock[ ::2, ::2], nImp/2, nEmb/2, RdmHl[ ::2, ::2], self.WfDecl.nElecA),
            FitComponent(EmbFock[1::2,1::2], nImp/2, nEmb/2, RdmHl[1::2,1::2], self.WfDecl.nElecB))
      pass

class FJobGroupResult(object):
   def __init__(self, FragmentResults):
      # sum up the total energy and total electron number
      TotalEnergy = 0.
      TotalElec = 0.
      for (Fragment, r) in FragmentResults:
         Factor = Fragment.GetTotalFactor()
         TotalEnergy += Factor * r.EmbEnergy
         TotalElec += Factor * r.EmbElec
      self.TotalElecEnergy = TotalEnergy
      self.TotalEnergy = None
      self.TotalElec = TotalElec
      self.FragmentResults = FragmentResults
      self.FullSystemFock = None
      self.FullSystemVloc = None

class FJobGroup:
   def __init__(self, FullSystem, Fragments):
      """represents the entire system and its subdivision into fragments."""
      assert(isinstance(FullSystem,FSystem))
      for F in Fragments:
         assert(isinstance(F,FFragmentDesc))
      self.FullSystem = FullSystem
      self.Fragments = Fragments
      self.CheckFragmentation()

   def CheckFragmentation(self):
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
         print "WARNING: Fragmenation covers no sites."
         Ok = False
      if ( MaxSite + 1 != self.FullSystem.nOrb ):
         print   "WARNING: Some orbitals of full system not present in fragmentation:"+\
               "\n         nOrb = %i  MaxSite = %i." % (self.FullSystem.nOrb, MaxSite)
         Ok = False
      for Site in range(MaxSite):
         if Site not in SiteOccurence:
            print "WARNING: Site %i not covered in fragmentation!" % Site
            Ok = False
         if abs(SiteOccurence[Site] - 1.) > 1e-8:
            print "WARNING: Site %i is covered %.5f times. Factor sum should be unity!" % (Site, CoveredSites[Site])
            Ok = False
      if not Ok:
         Msg = "\nERROR: There are problems with the fragmentation. Check and ignore."
         raise Exception(Msg)
         print Msg + "\n"

   def RunFragments(self, FullSystemHfResult, PrintEmb):
      """run HL electronic structure calculation on the individual
      fragments. return results of each of them."""
      FragmentResults = []
      for (iFragment,Fragment) in enumerate(self.Fragments):
         if PrintEmb:
            PrintHeading("EMBEDDED SUBSYSTEM: Id =%3i  F =%9.5f"\
               % (iFragment, Fragment.Factor), 0)
         sr = Fragment.Run(FullSystemHfResult, Print = PrintEmb)
         FragmentResults.append((Fragment, sr))
         if iFragment == STOP_AFTER_FRAGMENT:
            print "!! stopping after specified fragment."
            raise SystemExit
      return FJobGroupResult(FragmentResults)

   def FitCorrelationPotentials(self, JobGroupResult, FullSystemHfResult, PrintEmb):
      FragmentPotentials = []
      for (Fragment,r) in JobGroupResult.FragmentResults:
         vloc = Fragment.FitCorrelationPotential(r, FullSystemHfResult, Print=PrintEmb)
         FragmentPotentials.append(vloc)
      JobGroupResult.FragmentPotentials = FragmentPotentials

   def PrintFragmentResults(self, JobGroupResult, iMacroIter, dVsum, PrintEmb):
      """print summary of the fragment results."""
      jgr = JobGroupResult
      assert(isinstance(jgr,FJobGroupResult))

      PrintHeading("FRAGMENT SUMMARY & SYSTEM REASSEMBLY -- ITER. %i"\
         % (1+iMacroIter), int(PrintEmb))
      print ("  {:^5s} {:^13s} {:^14s}{:^12s}   {:>6s}").format(
         "SUB.SYS","FACTOR","ELEC.ENERGY", "ELECTRONS", "d[V]")
      for (iFragment, (Fragment,r)) in enumerate(jgr.FragmentResults):
         dV = "  --"
         if ( jgr.FragmentPotentials is not None ):
            dV = "%9.2e" % norm(jgr.FragmentPotentials[iFragment])
         print " %5i %13.5f %15.8f %11.6f %10s" % (iFragment,
            Fragment.Factor, r.EmbEnergy, r.EmbElec, dV)
      print
      print g.pResultFmt % ("Potential convergence", dVsum)
      print
      print g.pResultFmt % ("Number of electrons", jgr.TotalElec)
      print g.pResultFmt % ("Total electronic energy", jgr.TotalElecEnergy)
      print g.pResultFmt % ("Full.sys. core energy", self.FullSystem.CoreEnergy)
      print g.pResultFmt % ("Total energy", jgr.TotalEnergy)
      print

   def Run(self, VlocGuessLarge = None):
      PrintEmb = (g.IPRINT >= 2)
      FullCoreH_Real = 1. * self.FullSystem.CoreH
      if VlocGuessLarge is not None:
         VlocLarge = 1.*VlocGuessLarge
      else:
         VlocLarge = zeros(self.FullSystem.CoreH.shape)
      dc = [FDiisContext(4), None][0]
      ContinueDiis = False

      p = g.IPRINT == 0
      if p:
         PrintHeading("SELF-CONSISTENT EMBEDDING: CORR.POTENTIAL FIT", 0)
         Other = ""
         if hasattr(self.FullSystem.WfDecl,"SoBasis"):
            Other = " (spat irp := 0)"
         print "   THRDVC = %8.2e       VLOC_TYPE = %s/%s%s\n" % (THRDVC, g.VLOC_TYPE, g.VLOC_FIT_TYPE, Other)
         print ("   {:^5s} {:^15s} {:^14s}  {:>6s}  {:>6s}").format(
         "ITER.","TOT.ENERGY", "ELECTRONS", "d[V]", "DIIS")

      def RunVlocUpdate(VlocLarge):
         if PrintEmb: PrintHeading("FULL SYSTEM -- MEAN FIELD CALCULATION", 0)
         #if PrintEmb: PrintMatrix("VlocLarge before current it", VlocLarge)
         self.FullSystem.CoreH = FullCoreH_Real + VlocLarge
         FullSystemHfResult = self.FullSystem.RunHf(Print=PrintEmb)
         jgr = self.RunFragments(FullSystemHfResult, PrintEmb = (g.IPRINT >= 3))
         jgr.MeanFieldEnergy = FullSystemHfResult.Energy
         jgr.TotalEnergy = jgr.TotalElecEnergy + self.FullSystem.CoreEnergy
         self.FitCorrelationPotentials(jgr, FullSystemHfResult, PrintEmb)

         # now: add changed potentials to large system, and run again.
         dVlocLarge = zeros(VlocLarge.shape)
         if ( g.VLOC_TYPE == "Local" or g.VLOC_TYPE == "Diagonal" ):
            for (Fragment, VlocFull) in zip(self.Fragments, jgr.FragmentPotentials):
               nImp = len(Fragment.Sites)
               Vloc = Fragment.Factor * VlocFull[:nImp,:nImp]
               for Replica in Fragment.Replicas:
                  VlocR = mdot(Replica.Trafo, Vloc, Replica.Trafo.T)
                  E = list(enumerate(Replica.Sites))
                  for (i,iSite),(j,jSite) in it.product(E,E):
                     dVlocLarge[iSite,jSite] += VlocR[i,j]
         else:
            assert( g.VLOC_TYPE == "ImpAndEnv" )
            for (Fragment, VlocFull) in zip(self.Fragments, jgr.FragmentPotentials):
               assert(len(Fragment.Replicas) == 1)
               # ^- not impossible, but would be non-trivial to support.
               #    would require getting transformations of the entire basis,
               #    not just of the fragment sites.
               dVlocLarge += mdot(Fragment.EmbBasis, VlocFull, Fragment.EmbBasis.T)

         if hasattr(self.FullSystem.WfDecl,"SoBasis"): dVlocLarge = self.FullSystem.ProjectToIrrep0(dVlocLarge)
         return dVlocLarge, jgr, FullSystemHfResult

      Best_dVsum = 1e99
      try:
         for iMacroIter in range(MAXIT_VC):
            dVlocLarge, jgr, FullSystemHfResult = RunVlocUpdate(VlocLarge)
            dVsum = norm(dVlocLarge.flatten())
            if (iMacroIter == 0 ):
               InitialHfEnergy = FullSystemHfResult.Energy

            #if PrintEmb: print "dVlocLarge update:\nAlpha:\n%s\nBeta:\n%s\n" % ((dVlocLarge)[ ::2, ::2], (dVlocLarge)[1::2,1::2])
            if p:
               DiisText = " -  -"
               if dc: DiisText = "%2i %2i" % (dc.nDim, dc.iNext)
               print " %5i %15.8f %13.7f    %8.2e  %s" % (iMacroIter+1,
                  jgr.TotalEnergy, jgr.TotalElec, dVsum, DiisText)
            else:
               self.PrintFragmentResults(jgr, iMacroIter, dVsum, PrintEmb)

            if ( dVsum <= Best_dVsum ):
               jgr.FullSystemFock = FullSystemHfResult.GetFock()
               jgr.FullSystemVloc = VlocLarge
            jgr.dVc = dVsum
            if ( dVsum < THRDVC ):
               break

            Fock = FullSystemHfResult.GetFock() + dVlocLarge
            if dc is not None:
               SkipDiis = False
               #SkipDiis = dVsum > 1e-3 and iMacroIter < 5
               VlocLarge, dVlocLarge, Fock, c0 = dc.Apply(VlocLarge, dVlocLarge, Fock, Skip=SkipDiis)
               if not p: print " Vloc extrapolation: DIIS{%2i %2i  %8.2e}" % (dc.nDim, dc.iNext, c0)
            #self.FullSystem.FockGuess = Fock + dVlocLarge
            self.FullSystem.FockGuess = Fock
            VlocLarge += dVlocLarge
         if p and False:
            print
            PrintMatrix("Final correlation potential (diagonal)", diag(VlocLarge))
            PrintMatrix("Final mean field rdm (diagonal)", diag(FullSystemHfResult.GetRdm()))
         if p:
            print
            print g.pResultFmtAnnoted % ("Mean-field energy", InitialHfEnergy, "initial")
            print g.pResultFmtAnnoted % ("Mean-field energy", FullSystemHfResult.Energy, "final")
            print
            self.PrintFragmentResults(jgr, iMacroIter, dVsum, PrintEmb)
      finally:
         self.FullSystem.CoreH = self.FullSystem.CoreH_Unpatched

      #jgr.FullSystemFock = FullSystemHfResult.GetFock()
      #jgr.FullSystemVloc = VlocLarge
      if VlocLarge.shape[0] <= 12:
         print "Final VlocLarge:\n%s" % VlocLarge
      return jgr

# _____________________________________________________________________
#
#  Rest of the code deals with setup and execution of test calculations.
#



def MakeUnitRing(N):
   """makes a 3 x N position matrix representing a ring in which
   each two consecutive elements have unit distance from each other."""
   Pos = zeros((3,N))
   rH = .5/sqrt(sin((pi/N))**2)
   Phases = 2.*pi/N*arange(N)
   Pos[0,:] = rH * cos(Phases)
   Pos[1,:] = rH * sin(Phases)
   return Pos

def MakeUnitGrid2d(N,M):
   """Makes a 3 x N position matrix representing an equidistant N x M
   grid; lattice constant (i.e., nearest neighbor distance) is 1."""
   Pos = zeros((3,N*M))
   for i,j in it.product(range(N), range(M)):
      ij = i + N*j
      Pos[0,ij] = i - (N-1)/2.
      Pos[1,ij] = j - (M-1)/2.
   return Pos

def main():
   #OrbBasis = "STO-6G"
   OrbBasis = "MINAO"
   #OrbBasis = "cc-pVDZ"

   # stores how for each element the basis is defined in terms
   # of elementary functions.
   #    in []: core orbitals
   #    in (): polarization orbitals (usually uncontracted!)
   BasisDesc = MakeDefaultBasisDesc(OrbBasis)

   def PrintBasisDescForXyz(FileNameXyz):
      Xyz, AtomNames = ReadXyzFile(FileNameXyz,1./ToAng)
      Atoms = FAtomSet(Xyz,AtomNames)
      bdfa = BasisDesc.ForAtoms(Atoms)
      print "\n".join(map(str,bdfa))
      print "BasisLabels = [%s]" % ", ".join(['"%s"' % o for o in bdfa])
      print len(bdfa)
      raise SystemExit
   #PrintBasisDescForXyz("/home/cgk/calc/iron-complexes/fe_cn6.xyz")


   FitBasis = "def2-TZVPP-JKFIT"
   #FitBasis = "def2-QZVPP-MP2FI"
   BasisLibs = ["def2-nzvpp-jkfit.libmol", "minao.libmol", "emsl_sto6g.libmol", "emsl_cc-pVDZ.libmol"]
   PrintImportantLookingStuff = True
   #Method = "RHF"
   Method = "FCI"
   #Method = ("FCI0", "RHF")     # FCI as valence method, RHF as core method.
   #Method = ("MP2", "RHF")       # MP2 as valence method, RHF as core method.

   # these two may come in handy for symmetric systems:
   # nSitesPerFragment divides the entire site basis into consecutive
   # equal length fragments
   #nSitesPerFragment = 2
   nSitesPerFragment = 3
   #nSitesPerFragment = 5
   # ReplicateFragments tells the program to just explicitly calculate
   # one of them and copy the data to the other ones. Will obviously
   # give rubbish unless the fragments are equivlant (i.e.., doing that
   # in a ring is fine, doing it in a grid is not.).
   ReplicateFragments = False

   def MakeFragments(Atoms, Method, FullSystem):
      if ( type(Method) is tuple ):
         assert(len(Method) == 2)
         Method, MethodCore = Method
      else:
         MethodCore = Method

      Fragments = []

      if ( GeometryType == "xyz-phenazone" ):
         ToAtomList = lambda s: [int(o)-1 for o in s.split(',')]
         FragmentAtoms = map(ToAtomList, ["1,15,16,17", "8,19,20,21", "2,3,9,14,26", "6,7,18,4,5", "22,23,10,11", "12,13,24,25"])

         PrintHeading("CORE ENERGY DISTRIBUTION",0)
         print ("  {:^7s} {:^13s}   {:^14s}").format(
            "SUB.SYS","CORE ENERGY","ATOM LIST")
         CoreEnergySum = 0.
         for (iFragment,Fragment) in enumerate(FragmentAtoms):
            def FmtSite(i): return "%3i" % i
            FragmentCoreEnergy = Atoms.fSubsetCoreRepulsion(Fragment)
            print " {:5n} {:16.8f}    {:<s}".format(iFragment, FragmentCoreEnergy, " ".join(map(FmtSite,Fragment)))
            CoreEnergySum += FragmentCoreEnergy
         print
         print g.pResultFmt % ("Total core energy", CoreEnergySum)
         print
         if abs(CoreEnergySum - Atoms.fCoreRepulsion()) > 1e-10:
            raise Exception("Core energy fragmentation doesn't add up.")

         # make a list of basis offsets for all the atoms.
         BasisSizes = []
         BasisOffsets = [0]
         for Atom in Atoms.Elements:
            if Atom == "H": BasisSizes.append(1)
            elif Atom in ["C","N","O"]: BasisSizes.append(4 if g.ABSORB_HF_CORE else 5)
            else:
               assert(0)
            BasisOffsets.append(BasisOffsets[-1] + BasisSizes[-1])
         def ToSiteList(iAtoms):
            L = []
            for iAtom in iAtoms:
               L += list(range(BasisOffsets[iAtom], BasisOffsets[iAtom+1]))
            return L
         FragmentSites = map(ToSiteList, FragmentAtoms)
         for Sites in FragmentSites:
            Fragments.append(FFragment(Sites, 1., FullSystem, Method))
         return Fragments

      if ( GeometryType == "xyz" ):
         BaseFragmentLengths = []
         #print "x0:\n", Atoms.Orientations
         for (iAtom,Atom) in enumerate(Atoms.Elements):
            iElem = ElementNumbers[Atom]
            SplitCore = False and not g.ABSORB_HF_CORE  # no core left to split.
            #SplitCore = True
            # FIXME: do this properly, using BasisDesc.
            if ( Atom in ["H", "He"] ):
               if OrbBasis == "cc-pVDZ":
                  BaseFragmentLengths.append((5,eye(5),Method))  # 1s, 2s, 2px, 2py, 2pz
               else:
                  BaseFragmentLengths.append((1,eye(1),Method))  # 1s
            elif ( Atom in ["Li", "Be"] ):
               if g.ABSORB_HF_CORE:
                  BaseFragmentLengths.append((1,eye(1),Method))  # 2s -- with absorbed core.
               else:
                  if SplitCore:
                     BaseFragmentLengths.append((1,eye(1),MethodCore))  # 1s
                     BaseFragmentLengths.append((1,eye(1),Method))  # 2s
                  else:
                     BaseFragmentLengths.append((2,eye(2),Method))  # 1s, 2s
            elif ( Atom in "B C N O F Ne".split() ):
               if g.ABSORB_HF_CORE:
                  np = 4 # 2s+2p
               else:
                  if SplitCore:
                     BaseFragmentLengths.append((1,eye(1),MethodCore))  # 1s
                     np = 4 # 2s+2p
                  else:
                     np = 5 # 1s+2s+2p
               Trafo = eye(np)
               if ( Atoms.Orientations is not None ):
                  #print "x1:\n", Atoms.Orientations[:,:,iAtom].T
                  Trafo[(np-3):,(np-3):] = Atoms.Orientations[:,:,iAtom]
               BaseFragmentLengths.append((np,Trafo,Method))
               #assert(ReplicateFragments == False)
               ## ^- would require solid harmonic rotation!!
            elif ( Atom in "Br".split() ):
               assert(g.ABSORB_HF_CORE)
               BaseFragmentLengths.append((4,eye(4),Method))  # 1s
            elif ( Atom in "Fe".split() ):
               assert(g.ABSORB_HF_CORE)
               BaseFragmentLengths.append((5,eye(5),MethodCore))  # 4s + 5s + 4p
               BaseFragmentLengths.append((5,eye(5),Method))  # 3d
            else:
               assert(0)

         BaseFragmentSites = []
         BaseFragmentTrafos = []
         BaseFragmentMethods = []
         iSite = 0
         for (n,t,mthd) in BaseFragmentLengths:
            BaseFragmentSites.append(list(range(iSite,iSite+n)))
            BaseFragmentTrafos.append(t)
            BaseFragmentMethods.append(mthd)
            iSite += n

         ##for C-H fragments in benzene/MINAO:
         #for i in range(6):
            #Fragments.append(FFragment(list(range(5*i,5*(i+1))), 1., FullSystem, "FCI"))
            #Fragments.append(FFragment(list(range(5*i,5*(i+1))), 1., FullSystem, "CC(2)"))
         #for i in range(3):
            #Fragments.append(FFragment(list(range(10*i,10*(i+1))), 1., FullSystem, "CC(2)"))
         ##Fragments.append(FFragment(list(range(0,5)), 1., FullSystem, "CC(2)"))
         #Fragments.append(FFragment(list(range(0,10)), 1., FullSystem, "CC(2)"))
         ##Fragments.append(FFragment(list(range(5,10)), 1., FullSystem, "RHF"))
         #Fragments.append(FFragment(list(range(10,20)), 1., FullSystem, "RHF"))
         #Fragments.append(FFragment(list(range(20,30)), 1., FullSystem, "RHF"))
         ##Fragments.append(FFragment(list(range(10,20)), 1., FullSystem, "RHF"))
         ##Fragments.append(FFragment(list(range(20,30)), 1., FullSystem, "RHF"))
         #return Fragments

         # now: do other stuff with units of BaseFragmentSites instead
         # of individual basis functions.
         #    [...]

         nUniqueBaseFragments = [2,3][SplitCore]
         #nUniqueBaseFragments = [1,2][SplitCore]
         if ReplicateFragments:
            # calculate only a single unique fragment per type, but replicate it to other sites.
            for iUq in range(nUniqueBaseFragments):
               Clusters = range(iUq,len(BaseFragmentSites),nUniqueBaseFragments)
               Replicas = [FReplica(BaseFragmentSites[i], BaseFragmentTrafos[i]) for i in Clusters]
               Fragments.append(FFragment(BaseFragmentSites[iUq], 1.,
                  FullSystem, BaseFragmentMethods[iUq], Replicas=Replicas))
         else:
            # treat each fragment independently
            for (Sites, Method1) in zip(BaseFragmentSites,BaseFragmentMethods):
               Fragments.append(FFragment(Sites, 1., FullSystem, Method1))
         return Fragments

      N = FullSystem.nOrb
      assert(N % nSitesPerFragment == 0)
      Replicas = map(list,list(arange(N).reshape((N/nSitesPerFragment,nSitesPerFragment))))
      #Replicas = [ [0],[1,2,3], [4,5,6],[7,8,9] ]
      # ^- with diag fullrdm this almost correctly describes the habs-h10chain.
      #    problem seems to be the remaining 9 Hs, not the single H! Apparently
      #    chain is not correctly coupled to a doublet if not doing that.
      if ReplicateFragments:
         # calculate only a single fragment, but replicate it to other sites.
         Fragments.append(FFragment(Replicas[len(Replicas)/2], 1., FullSystem, Method, Replicas=Replicas))
      else:
         # treat each fragment independently
         for Sites in Replicas:
            Fragments.append(FFragment(Sites, 1., FullSystem, Method))
      return Fragments

   def RunAtomSet1(Atoms, Method, WfDecl, VlocGuess=None, FockGuess=None, GeomDesc = "", BasisDesc = None):
      if PrintImportantLookingStuff:
         PrintHeading("GEOMETRY AND ENVIRONMENT" + GeomDesc, 1)
         Scale,Units = 1., " [distances: a.u.]"
         Scale,Units = 0.5291772108, " [distances: angstrom]"
         print " Atomic configuration%s:\n" % Units
         print "   {:^5s}  {:^12s} {:^12s} {:^12s} {:^11s}".format(\
            "ELEM.","POS./X", "POS./Y", "POS./Z", "ORB.BASIS")
         for i in range(len(Atoms)):
            x,y,z = Scale * Atoms.Pos[:,i]
            print (" {:>5s}  {:12.6f} {:12.6f} {:12.6f} {:^16s}".format(\
               Atoms.Elements[i], x, y, z, OrbBasis)).rstrip()
         print
         print " %-34s" % "Fitting basis:" + FitBasis
         print
      FullSystem = MakeMolecularSystem(Atoms, WfDecl, OrbBasis, FitBasis, BasisLibs, BasisDesc)
      FullSystem.FockGuess = FockGuess

      if (g.WF_TYPE == "UHF" and FockGuess is None):
         print "\nWARNING: Hacking up special case Fock initial guess!"\
               "\n         strongly biasing initial occupation for a. .b a. .b a. .b ...\n"
         FockGuess = 1.*FullSystem.CoreH
         Bias = 1.
         #Bias = 0.01
         for i in xrange(FockGuess.shape[0]/4):
            FockGuess[4*i+0,4*i+0] -= Bias
            FockGuess[4*i+2,4*i+2] += Bias
            FockGuess[4*i+3,4*i+3] -= Bias
            FockGuess[4*i+1,4*i+1] += Bias
         FullSystem.FockGuess = FockGuess

         #if VlocGuess is None:
            #VlocGuess = .01*(FockGuess - FullSystem.CoreH)

      Fragments = MakeFragments(Atoms, Method, FullSystem)

      if PrintImportantLookingStuff and True:
         #print " Fragmentation:\n"
         PrintHeading("FRAGMENT AND REPLICA DEFINITION",0)
         print ("  {:^7s} {:^9s}{:^13s}   {:^14s}").format(
            "SUB.SYS","FACTOR","METHOD","SITE LIST")
         for (i,o) in enumerate(Fragments):
            def FmtSite(i): return "%3i" % i
            print " {:5n} {:11.5f} {:^13s}   {:<s}".format(i, o.Factor, o.Method, " ".join(map(FmtSite,o.Sites)))
            for Replica in o.Replicas:
               if ( tuple(Replica.Sites) == tuple(o.Sites) ):
                  continue
               print " {:5s} {:11s}  {:^10s}R -> {:<s}".format("","","", " ".join(map(FmtSite,Replica.Sites)))
               #print Replica.Trafo
         print
      JobGroup = FJobGroup(FullSystem, Fragments)
      jgr = JobGroup.Run(VlocGuess)
      return jgr

   def RunAtomSet(Atoms, Method, WfDecl, VlocGuess=None, FockGuess=None, GeomDesc = "", BasisDesc = None):
      iTries = 0
      while True:
         try:
            jgr = RunAtomSet1(Atoms, Method, WfDecl, VlocGuess, FockGuess, GeomDesc, BasisDesc)
         except ERhfConvergenceFail,e:
            iTries += 1
            nTries = 4
            print "\n\nWARNING: %s\n         Restarting current point (failed point: %i of %i).\n\n" % (str(e), iTries, nTries)
            #if ( iTries < nTries ):
               #continue
            #else:
            raise
         break
      if iTries != 0:
         print "\n\nWARNING: Current point was unstable. Needed %i tried to get it through!\n\n" % iTries
      return jgr

   def ScaleGeometry(OrigAtoms, Parameter, Args):
      """make geometry for symmetric stretching: multiply all coordinates by 'Parameter'"""
      return FAtomSet(Parameter * OrigAtoms.Pos, OrigAtoms.Elements)

   def MoveAtoms(OrigAtoms, Parameter, Args):
      """Args is a list of tuples (iAtom,vDir). Each iAtom will be moved by
      Parameter * vDir."""
      NewPos = 1. * OrigAtoms.Pos
      L = [o[0] for o in Args]
      assert(len(set(L)) == len(L)) # each one should occur only once.
      for (iAtom, vDir) in Args:
         assert(abs(norm(vDir) - 1.)<1e-8)
         assert(iAtom < len(OrigAtoms))
         NewPos[:,iAtom] += Parameter * vDir
      return FAtomSet(NewPos, OrigAtoms.Elements)

   def RunGeometrySet(UnscaledAtoms, WfDecl, DisplaceFnAndArgs, Parameters, Desc, PesFile):
      DisplaceFn, DisplaceArgs = DisplaceFnAndArgs
      VlocGuess = None
      FockGuess = None
      Lines = [Desc,("   " + 5*" {:^16s}").format("DISTANCE","TOT.ENERGY","ELECTRONS","MF.ENERGY", "  DIAG(VCORR)...")]
      for (iParam, Param) in enumerate(Parameters):
         try:
            Atoms = DisplaceFn(UnscaledAtoms, Param, DisplaceArgs)
            #g.VLOC_FIT_TYPE = "FullRdm"
            jgr = RunAtomSet(Atoms, Method, WfDecl, VlocGuess, FockGuess=FockGuess,
               GeomDesc = " - POINT %i of %i  [F =%8.5f]" % (iParam+1, len(Parameters), Param),
               BasisDesc = BasisDesc)
            VlocGuess = jgr.FullSystemVloc
            FockGuess = jgr.FullSystemFock
            #g.VLOC_FIT_TYPE = "ImpRdm"
            #jgr = RunAtomSet(Atoms, Method, WfDecl, VlocGuess, FockGuess=FockGuess,
               #GeomDesc = " - POINT %i of %i  [F =%8.5f]" % (iParam+1, len(Parameters), Param),
               #BasisDesc = BasisDesc)
            VlocDiag = diag(jgr.FullSystemVloc)
            if g.WF_TYPE == "UHF":
               v1 = 1.*VlocDiag
               VlocDiag[ ::2] = v1[::2] + v1[1::2]
               VlocDiag[1::2] = v1[::2] - v1[1::2]
            VlocDiag = " ".join(["%10.5f" % o for o in VlocDiag])
            Lines.append(" %16.8f %16.8f %16.8f %16.8f   %s" % (Param, jgr.TotalEnergy, jgr.TotalElec, jgr.MeanFieldEnergy, VlocDiag))
            PesText = "\n".join(Lines) + "\n"
            WriteFile(PesFile, PesText)
            print "*wrote PES to '%s'" % PesFile
            print
         except ERhfConvergenceFail,e:
            print "RHF convergence failed in final try. Skipping current point and resetting guess."
            raise
         except Exception,e:
            print "WARNING: Exception '%s' occurred." % e
            print "ignored!"
            raise
            pass
            #if VlocGuess is not None:
               #FockGuess -= VlocGuess
            #VlocGuess = None
            #FockGuess = None

      #g.VLOC_TYPE = "Local"
      #for (iParam, Param) in enumerate(Parameters[::-1]):
         #try:
            #Atoms = DisplaceFn(UnscaledAtoms, Param, DisplaceArgs)
            #jgr = RunAtomSet(Atoms, Method, WfDecl, VlocGuess, FockGuess=FockGuess,
               #GeomDesc = " - POINT %i of %i  [F =%8.5f]" % (iParam+1, len(Parameters), Param),
               #BasisDesc = BasisDesc)
            #VlocDiag = diag(jgr.FullSystemVloc)
            #if g.WF_TYPE == "UHF":
               #v1 = 1.*VlocDiag
               #VlocDiag[ ::2] = v1[::2] + v1[1::2]
               #VlocDiag[1::2] = v1[::2] - v1[1::2]
            #VlocDiag = " ".join(["%10.5f" % o for o in VlocDiag])
            #Lines.append(" %16.8f %16.8f %16.8f %16.8f   %s" % (Param, jgr.TotalEnergy, jgr.TotalElec, jgr.MeanFieldEnergy, VlocDiag))
            #VlocGuess = jgr.FullSystemVloc
            #FockGuess = jgr.FullSystemFock
            #PesText = "\n".join(Lines) + "\n"
            #WriteFile(PesFile, PesText)
            #print "*wrote PES to '%s'" % PesFile
            #print
         #except ERhfConvergenceFail,e:
            #print "RHF convergence failed in final try. Skipping current point and resetting guess."
            #raise
         #except Exception,e:
            #pass

      #print PesText
      pass

   rHH = 1.0 # given in a.u.
   GeometryType = "grid"
   #GeometryType = "ring"
   #GeometryType = "chain"
   #GeometryType = "xyz-phenazone"
   #GeometryType = "chain-habs"
   #GeometryType = "chain-h2abs"
   def MakeGeometry(GeometryType):
      Charge = 0
      SoBasis = None
      if GeometryType == "ring":
         nHyd = 10
         #nHyd = 6
         FileDesc = "-H-H-.x%i-ring" % nHyd
         Desc = "* --H--H-- (x%i ring); %s" % (nHyd, Method)
         Atoms = FAtomSet(rHH * MakeUnitRing(nHyd), nHyd * ["H"])
         rMin, rMax = 1.2, 6.
         DisplaceFnAndArgs = (ScaleGeometry, None)
      elif GeometryType == "grid":
         nHydX, nHydY = 4,3
         if ( nSitesPerFragment == 3 ):
            nHydX, nHydY = 3,4
         nHydX_, nHydY_ = max(nHydX,nHydY), min(nHydX, nHydY)
         FileDesc = "-H-H-.%ix%i-grid" % (nHydX_, nHydY_)
         Desc = "* --H--H-- (%ix%i grid); %s" % (nHydX_, nHydY_, Method)
         Atoms = FAtomSet(rHH * MakeUnitGrid2d(nHydX, nHydY), (nHydX*nHydY) * ["H"])
         rMin, rMax = 1.5, 6.0
         DisplaceFnAndArgs = (ScaleGeometry, None)

         if 1 and (nHydX, nHydY) == (3,4) or (nHydX, nHydY) == (4,3):
            # this may look like a very ugly hack, and to a certain degree it is.
            # Background: If things get complicated, resmin may run into trouble
            # and not be able to optimize its vloc hard enough; in that case
            # small symmetry contaminations can occur which confuse RHF
            # (it sees non-zero gradients, but cannot correct them because it can
            #  only really optimize in the symmetric subspace).
            # We work around that here by projecting dvloc and the Fock matrices
            # onto the totally symmetric spatial irrep, in which they are
            # supposed to lie.
            SoBasisX = array([[1,0,1],[0,1,0],[-1,0,1]]).T
            SoBasisY = array([[1,0,0,1],[0,1,1,0],[1,0,0,-1],[0,1,-1,0]]).T
            if (nHydX, nHydY) == (3,4):
               SoBasisX, SoBasisY = SoBasisY, SoBasisX
            nSo = SoBasisX.shape[1]*SoBasisY.shape[1]
            SoBasis = einsum("xi,yj->xyij",SoBasisX,SoBasisY).reshape(12,nSo)
            # normalize.
            Norms = diag(dot(SoBasis.T,SoBasis))**.5
            SoBasis = dot(SoBasis, diag(1./Norms))
            assert(allclose(dot(SoBasis.T,SoBasis) - eye(nSo),0.0))
            SoBasis = SoBasis[:,array([0,1,3,4,6,7,9,10,2,5,8,11])]
            #SoBasis = (SoBasis,array([4,4,2,2])) # <- basis, irrep lengths
            SoBasis = (SoBasis[:,:4],SoBasis[:,4:8],SoBasis[:,8:10],SoBasis[:,10:12])
            print "WARNING: restricting solution space to spatial irrep 0 orbitals in grid calculation."
            #raise SystemExit
      elif GeometryType == "chain":
         nHydX = 50
         assert(OrbBasis == "STO-6G")
         FileDesc = "-H-H-.x%i-chain" % (nHydX)
         Desc = "* --H--H-- (%i chain); %s" % (nHydX, Method)
         Xyz = rHH * MakeUnitGrid2d(nHydX, 1)
         Atoms = FAtomSet(rHH * Xyz, (nHydX*1) * ["H"])
         rMin, rMax = 1.2, 6.0
         DisplaceFnAndArgs = (ScaleGeometry, None)
         if 1:
            # see above, notes on grid.
            nSym = (nHydX+1)/2
            nAsym = nHydX/2
            assert(nSym+nAsym == nHydX)
            SoBasis = eye(nHydX)
            SoBasis[:,:nSym] = 2**(-.5) * (SoBasis[:,:nSym] + SoBasis[::-1,:nSym])
            SoBasis[:,nSym:] = 2**(-.5) * (SoBasis[:,nSym:] - SoBasis[::-1,nSym:])
            assert(allclose(dot(SoBasis.T,SoBasis) - eye(SoBasis.shape[0]),0.0))
            SoBasis = (SoBasis[:,:nSym],SoBasis[:,nSym:])
            print "WARNING: restricting solution space to spatial irrep 0 orbitals in chain calculation."
      elif GeometryType == "chain-habs":
         nHydX = 10
         FileDesc = "-H-H-.x%i-chain-habs" % (nHydX)
         Desc = "* --H--H-- (%i chain, H abstraction); %s" % (nHydX, Method)
         Xyz = 2.0 * MakeUnitGrid2d(nHydX, 1)
         # ^- 2.0 is about the equilibrium distance of the H10 chain.
         Atoms = FAtomSet(Xyz, (nHydX*1) * ["H"])
         rMin, rMax = -1.0,6.0
         Args = [(0, array([-1.,0.,0.]))]  # move first atom to the left.
         DisplaceFnAndArgs = (MoveAtoms, Args)
      elif GeometryType == "chain-h2abs":
         nHydX = 10
         FileDesc = "-H-H-.x%i-chain-h2abs" % (nHydX)
         Desc = "* --H--H-- (%i chain, H2 abstraction); %s" % (nHydX, Method)
         Xyz = 2.0 * MakeUnitGrid2d(nHydX, 1)
         # ^- 2.0 is about the equilibrium distance of the H10 chain.
         Atoms = FAtomSet(Xyz, (nHydX*1) * ["H"])
         rMin, rMax = -1.0,6.0
         Args = [(0, array([-1.,0.,0.])), (1, array([-1.,0.,0.]))]
         # ^- move first two atoms to the left.
         DisplaceFnAndArgs = (MoveAtoms, Args)
      elif GeometryType == "xyz-phenazone":
         FileNameXyz = "phenazone.xyz"; Charge=0; rMin, rMax= 1.,1.; Scale = 1./ToAng
         BaseName = path.splitext(path.basename(FileNameXyz))[0]
         DisplaceFnAndArgs = (ScaleGeometry, None)
         Xyz, AtomNames = ReadXyzFile(FileNameXyz,Scale)
         Atoms = FAtomSet(Xyz, AtomNames)
         Desc = '%s' % BaseName
         FileDesc = BaseName
      elif GeometryType == "xyz":
         #FileNameXyz = "/home/cgk/molecules/benzene.xyz"; rMin, rMax= 0.7, 3.0; Scale = 1./ToAng
         #FileNameXyz = "/home/cgk/molecules/li-ring-10.xyz"; rMin, rMax= 3.,20.; Scale = 1.
         #FileNameXyz = "/home/cgk/molecules/li-ring-4.xyz"; rMin, rMax= 3.,20.; Scale = 1.
         #FileNameXyz = "/home/cgk/molecules/li-ring-6.xyz"; rMin, rMax= 5.,10.; Scale = 1.
         #FileNameXyz = "/home/cgk/molecules/b-ring-6.xyz"; rMin, rMax= 2.,6.; Scale = 1.
         #FileNameXyz = "/home/cgk/molecules/h_pyramid_3.xyz"; rMin, rMax= 2.,6.; Scale = 1.
         #FileNameXyz = "/home/cgk/molecules/thomas/ch4.xyz"; rMin, rMax= 0.7, 3.0; Scale = 1./ToAng
         #FileNameXyz = "/home/cgk/calc/iron-complexes/fe_br4.xyz"; Charge=-2; rMin, rMax= 1.,1.; Scale = 1./ToAng
         FileNameXyz = "/home/cgk/calc/iron-complexes/fe_cn6.xyz"; Charge=-4; rMin, rMax= 1.,1.; Scale = 1./ToAng
         Xyz, AtomNames = ReadXyzFile(FileNameXyz,Scale)
         O = zeros((3,3,len(AtomNames)))
         for i in range(len(AtomNames)):
            eCenter = array([0.,0.,0.])
            eUp = array([0.,0.,1.])
            eFront = (Xyz[:,i] - eCenter); eFront /= norm(eFront)
            eRight = cross(eFront, eUp); eRight /= norm(eRight)
            eUp = cross(eRight, eFront); eUp /= norm(eUp)
            #O[:,:,i] = eye(3)
            O[:,0,i] = eRight
            O[:,1,i] = eFront
            O[:,2,i] = eUp
            #if ( AtomNames[i] == "C" ):
               #print "Atom pos: %s\nOrientation:\n%s\n" % (Xyz[:,i],O[:,:,i])

         DisplaceFnAndArgs = (ScaleGeometry, None)
         Atoms = FAtomSet(Xyz, AtomNames, Orientations=O)
         Desc = '%s' % BaseName
         FileDesc = BaseName
      else:
         raise Exception("GeometryType not recognized: '%s'" % GeometryType)

      if len(Atoms) != 50:
         assert(OrbBasis != "STO-6G")
      return Atoms, Desc, FileDesc, DisplaceFnAndArgs, rMin, rMax, Charge, SoBasis
   Atoms,Desc,FileDesc,DisplaceFnAndArgs,rMin,rMax,Charge,SoBasis = MakeGeometry(GeometryType)
   WfDecl = FWfDecl(nElec=Atoms.nElecNeutral()-Charge)
   if SoBasis is not None:
      WfDecl.SoBasis = SoBasis
   Desc += " // %s/%s // Core: %s" % (g.VLOC_TYPE, g.VLOC_FIT_TYPE, ["Included","AbsorbedViaHF"][g.ABSORB_HF_CORE])

   if ( nSitesPerFragment != 1 ):
      FileDesc += ".c%i" % nSitesPerFragment
   if not g.USE_INT2E:
      FileDesc += ".2eF"
   if g.WF_TYPE != "RHF":
      FileDesc += ".%s" % g.WF_TYPE

   if 1:
      Atoms.Pos *= 2.5
      RunAtomSet(Atoms, Method, WfDecl, BasisDesc=BasisDesc)
   else:
      #nDistances = 250
      nDistances = 50
      #nDistances = 100
      def MakeDistSet(Min, Max, N, Power):
         # this leads to points getting farer apart on the max side.
         return Min + (Max-Min) * arange(0.0, 1.0+1e-10, 1./N)**Power
      Direction = {"RHF":+1, "UHF":-1}[g.WF_TYPE]
      Power = 2.
      #Power = 1.
      RunGeometrySet(Atoms, WfDecl, DisplaceFnAndArgs, MakeDistSet(rMin, rMax, nDistances, Power)[::Direction], Desc,\
         "/home/cgk/dev/molecule_dmet/plots/pes.%s.txt" % FileDesc)

if __name__ == "__main__":
   iArg = 1
   while iArg < len(argv):
      Cmd = argv[iArg]
      if ( Cmd == "--iprint" ):
         iArg += 1
         g.IPRINT = int(argv[iArg])
      elif ( Cmd == "--fit" ):
         iArg += 1
         g.VLOC_FIT_TYPE = argv[iArg]
      else:
         raise Exception("Command line argument not recognized: %s" % Cmd)
      iArg += 1
   main()


