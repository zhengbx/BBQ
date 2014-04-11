import numpy as np
import scipy.linalg as la
from tempfile import mkdtemp
from os import remove, path
from commands import getoutput
from basicForNormal import *
import settings as dmetSet
from helpers import *

class FMfdResults(object):
   def __init__(self,Rdm,Mu,Gap):
      self.Gap = Gap
      self.Mu = Mu
      self.Rdm = Rdm

class FHlResults(object):
   def __init__(self, Energy, Energy2e, nElec, Rdm, cmd, output):
      self.Energy = Energy
      self.Energy2e = Energy2e
      self.nElec = nElec
      self.Rdm = Rdm
      self.cmd = cmd
      self.output = output

class NormalDmet(object):
   def __init__(self, OrbType, nElec=None, nElecA=None, nElecB=None, Ms2=None, ThrDeg=1.0e-8, ThrBathSvd = 1.0e-8):
      self.MfdThrDeg = ThrDeg
      self.MebThrBathSvd = ThrBathSvd
      self.OrbType = OrbType
      if nElecA is not None :
         self.MfdnElecA = nElecA
         self.MfdnElecB = nElecB
         self.MfdnElec = self.MfdnElecA + self.MfdnElecB
         self.MfdMs2 = self.MfdnElecA - self.MfdnElecB
      elif self.OrbType == "UHF" :
         raise Exception("number of alpha and beta electron must be specified")
      elif self.OrbType == "RHF" :
         assert (nElec%2==0)
         self.MfdnElec = nElec
         self.MfdMs2 = 0
         self.MfdnElecA = nElec/2
         self.MfdnElecB = nElec/2
      elif self.OrbType == "ROHF":
         raise Exception("restricted open shell hartree fock haven't implemented yet")


   def RunMfd(self, MfdHam, HamBlockDiag=False):
      if HamBlockDiag == True :
         if self.OrbType == "RHF":
            Ews, Orbs = Diag1eGQNHamiltonian(MfdHam)  
            Rdm, Mu, Gap = DealGQNOrbitals(Ews, Orbs, self.MfdnElecA, self.MfdThrDeg)  
            Rdm = 2*Rdm
         elif self.OrbType == "UHF":
            EwsA, OrbsA = Diag1eGQNHamiltonian(ExtractSpinComp(MfdHam,0))  
            RdmA, MuA, GapA = DealGQNOrbitals(EwsA, OrbsA, self.MfdnElecA, self.MfdThrDeg)  
            EwsB, OrbsB = Diag1eGQNHamiltonian(ExtractSpinComp(MfdHam,1))  
            RdmB, MuB, GapB = DealGQNOrbitals(EwsB, OrbsB, self.MfdnElecB, self.MfdThrDeg)  
            Rdm = CombineSpinComps(RdmA, RdmB)
            Mu = CombineSpinComps(MuA, MuB)
            Gap = CombineSpinComps(GapA, GapB)
      else :
         if self.OrbType == "RHF":
            Rdm, Mu, Gap = Diag1eHamiltonian(MfdHam, self.MfdnElec, self.MfdThrDeg) 
            nElecFull = np.trace(Rdm)
            #print "nElec is:" , nElecFull, int(nElecFull+0.1)
            #Rdm = 2*Rdm
            #nElecFull = np.trace(Rdm)
            #print "nElec is:" , nElecFull, int(nElecFull+0.1)
	    #raise SystemExit
         elif self.OrbType == "UHF":
            print MfdHam.shape
            RdmA, MuA, GapA = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecA, self.MfdThrDeg) 
            raise SystemExit
            RdmB, MuB, GapB = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecB, self.MfdThrDeg) 
            Rdm = CombineSpinComps(RdmA, RdmB)
            Mu = CombineSpinComps(MuA, MuB)
            Gap = CombineSpinComps(GapA, GapB)
      return FMfdResults(Rdm,Mu,Gap)


   def MakeEmbBasis(self,MfdRdm,ImpSites):
      if self.OrbType == "RHF":
         EmbBasis = MakeEmbeddingBasis(ImpSites, MfdRdm, self.MebThrBathSvd)
      elif self.OrbType == "UHF":
         EmbBasisA = MakeEmbeddingBasis(ImpSites[::2]/2, ExtractSpinComp(MfdRdm,0), self.MebThrBathSvd)
         EmbBasisB = MakeEmbeddingBasis(ImpSites[1::2]/2, ExtractSpinComp(MfdRdm,1), self.MebThrBathSvd)
         EmbBasis = CombineSpinComps(EmbBasisA, EmbBasisB)
      return EmbBasis

   def MakeEmbHam(self,MfdHam, MfdRdm, HAM, EmbBasis):
      EmbFock = ToEmb(MfdHam,EmbBasis)
      EmbRdm = ToEmb(MfdRdm,EmbBasis)
      #FIXME: need to change the 2el integrals here as well (for FCI input)
      #nElecFull = int(np.trace(MfdRdm))
      #print "nElec is:" , nElecFull
      return EmbFock, EmbFock, EmbRdm
         

   def ImpSolver(self,EmbFock,EmbRdm,HlMethod,nSysTrace=None,Int2e_=None,U=None):
      HlExecutable = {"Fci": dmetSet.FciExecutable,
                      "Dmrg": dmetSet.DmrgExecutable,
                      "CC": dmetSet.CCExecutable
                     }
      BasePath = mkdtemp(prefix=HlMethod, dir=dmetSet.TmpDir)
      HlInputName = path.join(BasePath, HlMethod+"INP")
      FileNameRdm = path.join(BasePath, HlMethod+"1RDM")
      nElec = int(np.trace(EmbRdm)+0.005)
      print "nElec is:" , nElec
      if U is not None :
         HlOption = OptionForHlMethod(HlMethod,FileNameRdm,BasePath,self.OrbType,nSys=EmbFock.shape[0], U=U, nSysTrace = nSysTrace)
         if self.OrbType == "RHF":
            WriteFciInput(HlInputName, EmbFock, nElec, U=1.0 )
         elif self.OrbType == "UHF":
            WriteFciInput(HlInputName, EmbFock, nElec, U=1.0, Uhf=True )
      else :
         assert(Int2e_ is not None)
         HlOption = OptionForHlMethod(HlMethod,FileNameRdm,BasePath,self.OrbType,nSys=EmbFock.shape[0], nSysTrace = nSysTrace)
         if self.OrbType == "RHF":
            Int2e = Int2e_
            WriteFciInput(HlInputName, EmbFock, nElec, Int2e=Int2e )
         elif self.OrbType == "UHF":
            Int2e = (Int2e_[0::2,0::2,0::2,0::2],
                     Int2e_[1::2,1::2,1::2,1::2],
                     Int2e_[0::2,0::2,1::2,1::2])
            WriteFciInput(HlInputName, EmbFock, nElec, Int2e=Int2e, Uhf=True )
      BaseCmd = HlExecutable[HlMethod]
      BaseCmd += HlOption
      Cmd = "%s '%s'" % (BaseCmd, HlInputName)
      print Cmd
      Output = getoutput(Cmd)
      print Output
      try:
         Energy = FindOutput(Output,"ENERGY")
         EpTrace = None
         if nSysTrace is not None:
            EpTrace = FindOutput(Output,"pTraceSys")
      except IndexError,e:
          print "Calculation failed: Output = ", Output
          print "Command was:"
      if ( self.OrbType == "RHF" ):
          Rdm1 = Read1Rdm(FileNameRdm)
      else :
          Rdm1a = Read1Rdm(FileNameRdm + ".A")
          Rdm1b = Read1Rdm(FileNameRdm + ".B")
          Rdm1 = CombineSpinComps(Rdm1a, Rdm1b)
      print "Rdm1 is:\n", Rdm1
      Cleanup(BasePath)
      return FHlResults(Energy, EpTrace, nElec, Rdm1, Cmd, Output)

   def GuessVcor(self, OrbType, dmet_inp, sc, U = 0., shift = 0.):
     # U could be hubbard U, in other cases, it's the upperlimit of Vcor entries
     # for some cases, though, U doesn't appear
     norbs = sc.nsites
     if OrbType == "UHF":
        norbs *= 2

     if dmet_inp.init_guess_type is None:
        Vcor = np.zeros((norbs, norbs))
        Vcor += np.eye(norbs) * shift
     elif dmet_inp.init_guess_type == 'MAN' and dmet_inp.init_guess is not None:
        Vcor = dmet_inp.init_guess
     elif dmet_inp.init_guess_type == "RAND":
        sites_with_V = []
        for frag in sc.fragments:
           if frag.get_emb_method() is not None: # none means do not do HL calculation
              sites_with_V += frag.sites.tolist()
        if OrbType == "UHF":
           sites_with_V = [2*x for x in sites_with_V] + [2*x+1 for x in sites_with_V]
        Vcor_small = np.random.rand(len(sites_with_V, sites_with_V)) * U/2
        Vcor_small += Vcor_small.T
        Vcor_small += np.eye(len(sites_with_V)) * shift
        Vcor[sites_with_V][:, sites_with_V] = Vcor_small
        return Vcor
     elif dmet_inp.init_guess_type == 'AF' and OrbType == "UHF":
        print "Warning: Automatic assignment of AF order may not be physical"
        assert(len(sc.fragments) == 1)
        Vcor_diag = np.array(norbs)
        Vcor_diag[::4] = U
        Vcor_diag[3::4] = U
        Vcor = np.diag(Vcor_diag)
        Vcor += np.eye(norbs) * shift
     else:
       raise Exception("Unknown or Unsupported Guess Type")

     return Vcor

if __name__ == '__main__':
    nElec = 10
    nOrb = 20
    Type = NormalDmet("RHF", nElec)
    MfdHam = np.random.random((nOrb,nOrb))
    MfdHam = MfdHam + MfdHam.T
    mfdResult = Type.RunMfd(MfdHam)
    EmbBasis = Type.MakeEmbBasis(mfdResult.Rdm,[0,4])
    EmbFock, EmbRdm = Type.MakeEmbHam(MfdHam, mfdResult.Rdm, EmbBasis) 
    hlResult = Type.ImpSolver(EmbFock,EmbRdm,"Fci",nSysTrace=2,U=5)
    #print hlResult.Energy
    #print hlResult.Energy2e




