import numpy as np
import scipy.linalg as la
from tempfile import mkdtemp
from os import remove, path
from commands import getoutput

from basicForNormal import *
import settings as dmetSet
from vcor_fit import FitVcorComponent
from helpers import *

class FMfdResults(object):
    def __init__(self,Rdm,Mu,Gap, Energy):
        self.Gap = Gap
        self.Mu = Mu
        self.Rdm = Rdm
        self.Energy = Energy       

class FHlResults(object):
    def __init__(self, Energy, Energy2e, nEmbElec, EmbEnergy, Rdm, cmd, output):
        self.Energy = Energy
        self.Energy2e = Energy2e
        self.EmbEnergy = EmbEnergy
        self.nEmbElec = nEmbElec
        self.Rdm = Rdm
        self.cmd = cmd
        self.output = output

class NormalDmet(object):
    def __init__(self, inp_mfd, inp_fitting, OrbType, nElec=None, nElecA=None, nElecB=None, Ms2=None):
        self.MfdThrDeg = inp_mfd.thr_deg
        self.MebThrBathSvd = inp_fitting.thr_svd_bath
        self.MfdAlgo = inp_mfd.mfd_algo
        self.MaxIter = inp_mfd.max_iter
        self.Thr_dE = inp_mfd.conv_thr_dE
        self.Thr_dOrb = inp_mfd.conv_thr_dOrb
        self.diis_thr = inp_mfd.diis_thr
        self.diis_start = inp_mfd.diis_start
        self.diis_dim = inp_mfd.diis_dim
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
            raise Exception("restricted open-shell hartree-fock not implemented yet")

    def RunMfd(self, MfdHam, Int2e, Lattice, HamBlockDiag=False):
        if HamBlockDiag == True :
            if self.MfdAlgo == 'diag1el':
                if self.OrbType == "RHF":
                    Ews, Orbs = Diag1eGQNHamiltonian(MfdHam)  
                    Rdm, Mu, Gap = DealGQNOrbitals(Ews, Orbs, self.MfdnElecA, self.MfdThrDeg)  
                    Rdm = 2*Rdm
                elif self.OrbType == "UHF":
                    EwsA, OrbsA = Diag1eGQNHamiltonian(ExtractSpinCompK(MfdHam,0))
                    RdmA, MuA, GapA = DealGQNOrbitals(EwsA, OrbsA, self.MfdnElecA, self.MfdThrDeg)  
                    EwsB, OrbsB = Diag1eGQNHamiltonian(ExtractSpinCompK(MfdHam,1))  
                    RdmB, MuB, GapB = DealGQNOrbitals(EwsB, OrbsB, self.MfdnElecB, self.MfdThrDeg)  
                    Rdm = CombineSpinCompsK(RdmA, RdmB)
                    Mu = np.array([MuA, MuB])
                    Gap = np.array([GapA, GapB])
                Energy = np.dot(Rdm.flatten(),MfdHam.flatten())
            elif self.MfdAlgo == 'hf':
                Rdm, Mu, Gap, Energy = RunHf_kspace(MfdHam, Int2e, self.MfdnElecA, self.MfdnElecB, 
                                           self.OrbType, Lattice, self.MfdThrDeg, self.MaxIter, 
                                           self.Thr_dE, self.Thr_dOrb, 
                                           self.diis_start, self.diis_thr, self.diis_dim)  
        else :
            if self.MfdAlgo == 'diag1el':
                if self.OrbType == "RHF":
                    Rdm, Mu, Gap = Diag1eHamiltonian(MfdHam, self.MfdnElec, self.MfdThrDeg) 
                    Rdm = 2*Rdm
                elif self.OrbType == "UHF":
                    RdmA, MuA, GapA = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecA, self.MfdThrDeg) 
                    RdmB, MuB, GapB = Diag1eHamiltonian(ExtractSpinComp(MfdHam,0), self.MfdnElecB, self.MfdThrDeg) 
                    Rdm = CombineSpinComps(RdmA, RdmB)
                    Mu = np.array([MuA, MuB])
                    Gap = np.array([GapA, GapB])
                    Energy = np.dot(Rdm.flatten(),MfdHam.flatten())
            elif self.MfdAlgo == 'hf':                   
                Rdm, Mu, Gap, Energy = RunHf_xspace(MfdHam, Int2e, self.MfdnElecA, self.MfdnElecB, 
                                            self.OrbType, Lattice, self.MfdThrDeg, self.MaxIter, 
                                            self.Thr_dE, self.Thr_dOrb, 
                                            self.diis_start, self.diis_thr, self.diis_dim)  

        return FMfdResults(Rdm,Mu,Gap, Energy)

    def MakeEmbBasis(self,MfdRdm,ImpSites):
        if self.OrbType == "RHF":
            EmbBasis = MakeEmbeddingBasis(ImpSites, MfdRdm, self.MebThrBathSvd)
        elif self.OrbType == "UHF":
            EmbBasisA = MakeEmbeddingBasis(ImpSites, np.around(ExtractSpinComp(MfdRdm,0),decimals=8), self.MebThrBathSvd)
            EmbBasisB = MakeEmbeddingBasis(ImpSites, np.around(ExtractSpinComp(MfdRdm,1),decimals=8), self.MebThrBathSvd)
            EmbBasis = CombineSpinComps(EmbBasisA, EmbBasisB)
        return EmbBasis
 
    def MakeEmbHam(self, Lattice, CoreH_aug, CoreH, MfdRdm, Int2e_, EmbBasis, fragsites):
        EmbCoreH_aug = ToEmb(Lattice, CoreH_aug, EmbBasis)
        EmbCoreH_ = ToEmb(Lattice, CoreH, EmbBasis)
        EmbCoreH_alpha = np.zeros_like(EmbCoreH_aug[ ::2, ::2])
        EmbCoreH_beta = np.zeros_like(EmbCoreH_aug[1::2,1::2])
        EmbCoreH_alpha += EmbCoreH_aug[ ::2, ::2]       
        EmbCoreH_beta += EmbCoreH_aug[1::2,1::2]
        for i in fragsites:
            for j in fragsites:
                 EmbCoreH_alpha[i,j] = EmbCoreH_[ ::2, ::2][i,j]
                 EmbCoreH_beta[i,j] = EmbCoreH_[1::2,1::2][i,j]
        EmbCoreH = np.zeros(EmbCoreH_aug.shape)
        EmbCoreH[ ::2, ::2] += EmbCoreH_alpha
        EmbCoreH[1::2,1::2] += EmbCoreH_beta
        
        EmbRdm = ToEmb(Lattice,  MfdRdm, EmbBasis)
        EmbDim = EmbBasis.shape[1]
        if Int2e_[0] is None:
            assert(Int2e_[1] is not None)
            U = Int2e_[1]
            Int2e_4ix = None
        else:
            assert(Int2e_[0] is not None)
            U = None
            Int2e_ = np.zeros(EmbDim)
            Int2e_4ix = np.zeros((EmbDim,EmbDim,EmbDim,EmbDim))
            # FIXME: insert trafo for true 2el interaction
        return EmbCoreH, EmbCoreH_aug, (Int2e_4ix, U), EmbRdm
          
 
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

    def ImpSolver(self, EmbCoreH, EmbRdm, HlMethod, n2eOrb, Int2e):
        Int2e_= Int2e[0]
        U_ = Int2e[1]
        BasePath = mkdtemp(prefix=HlMethod, dir=dmetSet.TmpDir)
        HlInputName = path.join(BasePath, HlMethod+"INP")
        nElec = int(np.trace(EmbRdm)+0.005)
        if self.OrbType == "UHF":
            Ms2_ = int(np.trace(EmbRdm[::2,::2])+0.005)-int(np.trace(EmbRdm[1::2,1::2])+0.005)
        if U_ is not None :
            assert(U_ is not None)
            if self.OrbType == "RHF":
                WriteDumpFile(HlInputName, EmbCoreH, nElec, n2eOrb, U=U_ )
            elif self.OrbType == "UHF":
                n2eOrb_ = [None] * 2 *len(n2eOrb)
                n2eOrb_[:len(n2eOrb)] = n2eOrb 
                n2eOrb_[len(n2eOrb):] = n2eOrb + len(n2eOrb)
                n2eOrb = n2eOrb_
                WriteDumpFile(HlInputName, EmbCoreH, nElec, n2eOrb, U=U_, Uhf=True, Ms2 = Ms2_ )
        else :
            assert(Int2e_ is not None)
            if self.OrbType == "RHF":
                WriteDumpFile(HlInputName, EmbCoreH, nElec, n2eOrb, Int2e=Int2e_)
            elif self.OrbType == "UHF":
                Int2e = (Int2e_[0::2,0::2,0::2,0::2],
                         Int2e_[1::2,1::2,1::2,1::2],
                         Int2e_[0::2,0::2,1::2,1::2])
                WriteDumpFile(HlInputName, EmbCoreH, nElec, n2eOrb, Int2e=Int2e_, Uhf=True, Ms2 = Ms2_)
        if HlMethod == "Fci":
            Energy, EpTrace, Rdm1, Cmd, Output = FindFciSolution(BasePath,np.shape(EmbCoreH)[0],len(n2eOrb),IntType=self.OrbType)
        elif HlMethod == "Dmrg":
            dmrgConf = WriteDmrgConfigFile(BasePath,np.shape(EmbCoreH)[0],nElec)
            Energy, EpTrace, Rdm1, Cmd, Output = FindDmrgSolution(BasePath,dmrgConf)
        else :
            raise Exception("You have to implement the interface for other solver")
        if n2eOrb is not None:
            nEmbElec = np.trace(Rdm1[:,n2eOrb][n2eOrb,:])
            # FIXME: insert EmbEnergy calc
            EmbEnergy = None
        return FHlResults(Energy, EpTrace, nEmbElec, EmbEnergy, Rdm1, Cmd, Output)
         
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
            Vcor_diag = np.zeros((norbs))
            Vcor_diag[::4] = U
            Vcor_diag[3::4] = U
            Vcor = np.diag(Vcor_diag)
            Vcor += np.eye(norbs) * shift
        else:
            raise Exception("Unknown or Unsupported Guess Type")

        return Vcor

    def FitCorrelationPotential(self, inp_fit, fragsites, nElec, EmbMfdHam, RdmHl):
        """fit a matrix nEmb x nEmb, referring to an operator connecting  
        all embedded sites to each other, to more closely approach the 
        correlated calculation result in with the mean-field method.
        Returns nEmb x nEmb matrix (always), even if due to the fitting
        criterion used the actual matrix is only a subset of this."""

        nImp = len(fragsites)
        # ^ FIXME: allow FitVcorComponent to specify the fragsites, not assuming them to be
        # the first sites in RmdHl
        VcorType = inp_fit.v_fit_domain
        VcorFitType = inp_fit.vfit_method
        if ( self.OrbType == "RHF" ):
            assert(nElec%2 == 0)
            return FitVcorComponent(EmbMfdHam, nImp, (1./2.)*RdmHl, VcorType, VcorFitType)
        else:
            assert( self.OrbType == "UHF")
            return CombineSpinComps(
                FitVcorComponent(EmbMfdHam[ ::2, ::2], nImp, RdmHl[ ::2, ::2], 
                                VcorType, VcorFitType),
                FitVcorComponent(EmbMfdHam[1::2,1::2], nImp, RdmHl[1::2,1::2], 
                                VcorType, VcorFitType))
        pass

    def VcorUpdate(self, FragmentPotentials, Fragments, shape):
        if self.OrbType == "RHF":
            for ( iFragment, Fragment) in enumerate(Fragments):
                VcorFull = FragmentPotentials[iFragment]
                fragsites = Fragment.get_sites()
                Vcor = Fragment.factor * VcorFull[fragsites,fragsites]
                for Replica in Fragment.get_replicas():
                    VcorR = mdot(Replica[1], Vcor, Replica[1].T)
                    E = list(enumerate(Replica[0]))
                    for (i,iSite),(j,jSite) in it.product(E,E):
                        dVcorLarge[iSite,jSite] += VcorR[i,j]
        elif self.OrbType == "UHF":
            dVcorLarge = np.zeros(shape)
            dVcorLarge_alpha = np.zeros(dVcorLarge[ ::2, ::2].shape)
            dVcorLarge_beta = np.zeros(dVcorLarge[1::2,1::2].shape)
            for ( iFragment, Fragment) in enumerate(Fragments):
                VcorFull = FragmentPotentials[iFragment]
                fragsites = Fragment.get_sites()
                Vcor = np.zeros((2*len(fragsites),2*len(fragsites)))
                Vcor[ ::2, ::2] = Fragment.factor * VcorFull[ ::2, ::2][fragsites,:][:,fragsites]
                Vcor[1::2,1::2] = Fragment.factor * VcorFull[1::2,1::2][fragsites,:][:,fragsites]
                for Replica in Fragment.get_replicas():
                    VcorR = mdot(Replica[1], Vcor[ ::2, ::2], Replica[1].T)
                    E = list(enumerate(Replica[0]))
                    for (i,iSite),(j,jSite) in it.product(E,E):
                        dVcorLarge_alpha[iSite,jSite] += VcorR[i,j]
                    VcorR = mdot(Replica[1], Vcor[1::2,1::2], Replica[1].T)
                    E = list(enumerate(Replica[0]))
                    for (i,iSite),(j,jSite) in it.product(E,E):
                        dVcorLarge_beta[iSite,jSite] += VcorR[i,j]
            dVcorLarge[ ::2, ::2] += dVcorLarge_alpha
            dVcorLarge[1::2,1::2] += dVcorLarge_beta
        return dVcorLarge            

if __name__ == '__main__':
    nElec = 10
    nOrb = 20
    Type = NormalDmet("RHF", nElec)
    MfdHam = np.random.random((nOrb,nOrb))
    MfdHam = MfdHam + MfdHam.T
    mfdResult = Type.RunMfd(MfdHam)
    EmbBasis = Type.MakeEmbBasis(mfdResult.Rdm,[0,4])
    EmbCoreH, EmbRdm = Type.MakeEmbHam(MfdHam, mfdResult.Rdm, EmbBasis) 
    nImp = 2
    #hlResult = Type.ImpSolver0(EmbCoreH,EmbRdm,"Fci",n2eOrb=nImp,U=1.0)
    hlResult = Type.ImpSolver(EmbCoreH,EmbRdm,"Fci",n2eOrb=nImp,U_=1.0)
    hlResult = Type.ImpSolver(EmbCoreH,EmbRdm,"Dmrg",n2eOrb=nImp,U_=1.0)
    #print hlResult.Energy
    #print hlResult.Energy2e
    




