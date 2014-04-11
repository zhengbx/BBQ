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
                 RdmA, MuA, GapA = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecA, self.MfdThrDeg) 
                 RdmB, MuB, GapB = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecB, self.MfdThrDeg) 
                 Rdm = CombineSpinComps(RdmA, RdmB)
                 Mu = CombineSpinComps(MuA, MuB)
                 Gap = CombineSpinComps(GapA, GapB)
        return FMfdResults(Rdm,Mu,Gap)

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
                RdmA, MuA, GapA = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecA, self.MfdThrDeg) 
                RdmB, MuB, GapB = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecB, self.MfdThrDeg) 
                Rdm = CombineSpinComps(RdmA, RdmB)
                Mu = CombineSpinComps(MuA, MuB)
                Gap = CombineSpinComps(GapA, GapB)
        return FMfdResults(Rdm,Mu,Gap)

    def MakeEmbHam(self,MfdHam, MfdRdm, HAM, EmbBasis):
        EmbFock = ToEmb(MfdHam,EmbBasis)
        EmbRdm = ToEmb(MfdRdm,EmbBasis)
        #FIXME: need to change the 2el integrals here as well (for FCI input)
        #nElecFull = int(np.trace(MfdRdm))
        #print "nElec is:" , nElecFull
        return EmbFock, EmbFock, EmbRdm

    def ImpSolver0(self,EmbFock,EmbRdm,HlMethod,n2eOrb,nSysTrace=None,Int2e_=None,U=None):
        HlExecutable = {"Fci": dmetSet.FciExecutable,
                        "Dmrg": dmetSet.DmrgExecutable,
                        "CC": dmetSet.CCExecutable
                        }
        BasePath = mkdtemp(prefix=HlMethod, dir=dmetSet.TmpDir)
        HlInputName = path.join(BasePath, HlMethod+"INP")
        FileNameRdm = path.join(BasePath, HlMethod+"1RDM")
        nElec = int(np.trace(EmbRdm)+0.005)
        if self.OrbType == "UHF":
            Ms2 = int(np.trace(EmbRdm[::2,::2])+0.005)-int(np.trace(EmbRdm[1::2,1::2])+0.005)
        print "nElec is:" , nElec
        if U is not None :
            HlOption = OptionForHlMethod(HlMethod,FileNameRdm,BasePath,self.OrbType,nSys=EmbFock.shape[0], U=U, nSysTrace = nSysTrace)
            if self.OrbType == "RHF":
                WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, U=U )
            elif self.OrbType == "UHF":
               WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, U=U, Uhf=True, Ms2 = Ms2 )
        else :
            assert(Int2e_ is not None)
            HlOption = OptionForHlMethod(HlMethod,FileNameRdm,BasePath,self.OrbType,nSys=EmbFock.shape[0], nSysTrace = nSysTrace)
            if self.OrbType == "RHF":
                Int2e = Int2e_
                WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, Int2e=Int2e )
            elif self.OrbType == "UHF":
                Int2e = (Int2e_[0::2,0::2,0::2,0::2],
                         Int2e_[1::2,1::2,1::2,1::2],
                         Int2e_[0::2,0::2,1::2,1::2])
            WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, Int2e=Int2e, Uhf=True, Ms2 = Ms2)
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
        print "fci Rdm1 is:\n", Rdm1
        Cleanup(BasePath)
        return FHlResults(Energy, EpTrace, nElec, Rdm1, Cmd, Output)

    def ImpSolver(self,EmbFock,EmbRdm,HlMethod,n2eOrb,Int2e_=None,U_=None):
        BasePath = mkdtemp(prefix=HlMethod, dir=dmetSet.TmpDir)
        HlInputName = path.join(BasePath, HlMethod+"INP")
        nElec = int(np.trace(EmbRdm)+0.005)
        if self.OrbType == "UHF":
           Ms2_ = int(np.trace(EmbRdm[::2,::2])+0.005)-int(np.trace(EmbRdm[1::2,1::2])+0.005)
        if U_ is not None :
           assert(U_ is not None)
           if self.OrbType == "RHF":
               WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, U=U_ )
           elif self.OrbType == "UHF":
               WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, U=U_, Uhf=True, Ms2 = Ms2_ )
        else :
            assert(Int2e_ is not None)
            if self.OrbType == "RHF":
                WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, Int2e=Int2e_)
            elif self.OrbType == "UHF":
                Int2e = (Int2e_[0::2,0::2,0::2,0::2],
                         Int2e_[1::2,1::2,1::2,1::2],
                         Int2e_[0::2,0::2,1::2,1::2])
            WriteDumpFile(HlInputName, EmbFock, nElec, n2eOrb, Int2e=Int2e_, Uhf=True, Ms2 = Ms2_)
        if HlMethod == "Fci":
            Energy, EpTrace, Rdm1, Cmd, Output = FindFciSolution(BasePath,np.shape(EmbFock)[0],n2eOrb)
        elif HlMethod == "Dmrg":
            dmrgConf = WriteDmrgConfigFile(BasePath,np.shape(EmbFock)[0],nElec)
            Energy, EpTrace, Rdm1, Cmd, Output = FindDmrgSolution(BasePath,dmrgConf)
        else :
            raise Exception("You have to implement the interface for other solver")
        
        return FHlResults(Energy, EpTrace, nElec, Rdm1, Cmd, Output)
          


if __name__ == '__main__':
    nElec = 10
    nOrb = 20
    Type = NormalDmet("RHF", nElec)
    MfdHam = np.random.random((nOrb,nOrb))
    MfdHam = MfdHam + MfdHam.T
    mfdResult = Type.RunMfd(MfdHam)
    EmbBasis = Type.MakeEmbBasis(mfdResult.Rdm,[0,4])
    EmbFock, EmbRdm = Type.MakeEmbHam(MfdHam, mfdResult.Rdm, EmbBasis) 
    nImp = 2
    #hlResult = Type.ImpSolver0(EmbFock,EmbRdm,"Fci",n2eOrb=nImp,U=1.0)
    hlResult = Type.ImpSolver(EmbFock,EmbRdm,"Fci",n2eOrb=nImp,U_=1.0)
    hlResult = Type.ImpSolver(EmbFock,EmbRdm,"Dmrg",n2eOrb=nImp,U_=1.0)
    #print hlResult.Energy
    #print hlResult.Energy2e
    




