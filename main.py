# Draft for the main driver routine for joint dmet code.
# It needs the classes: HAM, GEOM, TYPE scetched out below.
# Here, TYPE is variable that is assigned any of the allowed
# TYPE-classes: Normal, BCS, Eph.
# No print-out has been inserted yet.
# It is not yet connected with the input module.

import numpy as np
import itertools as it

from helpers import mdot, ExtractSpinComp, CombineSpinComps
from vcor_fit import FitVcorComponent
from results import FDmetResult
from diis import FDiisContext



# setup hamiltonian                                                                                
# provides: CoreH, CoreEnergy

# setup geometry
# provides: _EnergyFactor, Fragments, nImp,OrbType, ORB_OCC, nElec
# Fragment needs: Sites, Factor, Replicas, FullSystem, Method (class from Geralds code) 

# setup type classes
# provides: 
#    - MakeMfdHam: MfdHam 
#    - RunMfd: MfdResult (MeanFieldEnergy, MfdHam, CoreEnergy (class from Geralds code))
#    - MakeEmbedBasis: nEmb
#    - MakeEmbedHam:
#    - ImpSolver: Results(=FragmentResults)

# setup init guess (StartingGuess) 





def main(StartingGuess, GEOM, HAM, TYPE, INP):

    def RunFragments(Fragments, Rdm, MfdHam, HAM):
        FragmentResults = []
        for (iFragment,Fragment) in enumerate(Fragments):
            EmbBasis = TYPE.MakeEmbBasis(Rdm, Fragment)
            EmbHam, EmbMfdHam = TYPE.MakeEmbHam(EmbBasis, MfdHam, HAM, Fragment)
            Results = TYPE.ImpSolver(EmbHam, EmbMfdHam,Fragment)
            FragmentResults.append((Fragment, Results))
        ClusterFactor = GEOM._EnergyFactor()
        # ^- for normalization with super-cell cluster size. Purely cosmetic.
        return FDmetResults(FragmentResults, ClusterFactor)


    def FitCorrelationPotential(GEOM, TYPE, Fragment, MfdHam, VcorType, VcorFitType):
      """fit a matrix nEmb x nEmb, referring to an operator connecting  
      all embedded sites to each other, to more closely approach the 
      correlated calculation result in with the mean-field method.
      Returns nEmb x nEmb matrix (always), even if due to the fitting
      criterion used the actual matrix is only a subset of this."""
 
      MfdResult = TYPE.RunMfd(MfdHam, None)
      Rdm = MfdResult.Rdm
      EmbBasis = TYPE.MakeEmbBasis(Rdm, Fragment)
      EmbHam, EmbMfdHam = TYPE.MakeEmbHam(EmbBasis, MfdHam, Ham,Fragment)
      HlResults = TYPE.ImpSolver(EmbHam, EmbMfdHam,Fragment)
      RdmHl = HlResults.Rdm
      
      nImp = GEOM.nImp
      nEmb = TYPE.MakeEmbBasis.nEmb
      OrbType = GEOM.OrbType
      ORB_OCC = GEOM.OrbOcc()
      nElec = GEOM.nElec

      if ( OrbType == "RHF" ):
          assert(nElec%2 == 0)
          return FitVcorComponent(EmbMfdHam, nImp, (1./ORB_OCC)*RdmHl, VcorType, VcorFitType)
      else:
          assert(OrbType == "UHF")
          return CombineSpinComps(
              FitVcorComponent(EmbMfdHam[ ::2, ::2], nImp/2, RdmHl[ ::2, ::2], 
                               VcorType, VcorFitType),
              FitVcorComponent(EmbMfdHam[1::2,1::2], nImp/2, RdmHl[1::2,1::2], 
                               VcorType, VcorFitType))
      pass

    def FitCorrelationPotentials(GEOM, TYPE, MfdHam, VcorType, VcorFitType):
        FragmentPotentials = []
        for Fragment in Fragments:
            vloc = FitCorrelationPotential(GEOM, TYPE, Fragment, MfdHam, VcorType, VcorFitType)
            FragmentPotentials.append(vloc)
        return FragmentPotentials
    

    def RunVcorUpdate(VcorLarge, MfdHam, HAM):
        MfdResult = TYPE.RunMfd(MfdHam, VcorLarge)
        Rdm = MfdResult.Rdm
        FragRes = RunFragments(Fragments, Rdm, MfdHam, HAM)
        FragRes.MeanFieldEnergy = MfdResult.MeanFieldEnergy
        FragRes.TotalEnergy = FragRes.TotalElecEnergy + MfdResult.CoreEnergy
        FragRes.FragmentPotentials = FitCorrelationPotentials(GEOM, TYPE, MfdHam, VcorType, VcorFitType)
        dVcorLarge = np.zeros(VcorLarge.shape)
        assert(VcorType == "Local" or VcorType == "Diagonal")
        
        for (Fragment, VcorFull) in zip(GEOM.Fragments, FragRes.FragmentPotentials):
            nImp = len(Fragment.Sites)
            Vcor = Fragment.Factor * VcorFull[:nImp,:nImp]
            for Replica in Fragment.Replicas:
                VcorR = mdot(Replica.Trafo, Vcor, Replica.Trafo.T)
                E = list(enumerate(Replica.Sites))
                for (i,iSite),(j,jSite) in it.product(E,E):
                    dVcorLarge[iSite,jSite] += VcorR[i,j]
                    
        return dVcorLarge, FragRes, MfdResult


    #input parameters
    DmetMaxIt = INP.DMET.MaxIt
    ThrdVcor = INP.DMET.ThrdVcor
    DiisThr = INP.DMET.DiisThr
    DiisStart = INP.DMET.DiisStart
    DiisDim = INP.DMET.DiisDim
    VcorType = INP.FIT.VcorType
    VcorFitType = INP.FIT.VcorFitType


    MfdHam = np.zeros_like(TYPE.MfdHam)
    MfdHam += TYPE.MfdHam

    Fragments = GEOM.Fragments

    dc = FDiisContext(DiisDim)

    if 'vcor' in StartingGuess:
        VcorLarge = 1.*StartingGuess['vcor']
    else:
        #replace with data from GEOM
        VcorLarge = np.zeros_like(Ham[0,:,:])
          

    for iMacroIt in range(DmetMaxIt):
        dVcorLarge, FragRes, FSMfdResult = RunVcorUpdate(VcorLarge, MfdHam, HAM)
        dVsum = np.dot(dVcorLarge.flatten(),dVcorLarge.flatten())
        if (iMacroIt == 0 ):
            InitialHfEnergy = FSMfdResult.MeanFieldEnergy
        FragRes.dVc = dVsum
        if ( dVsum < ThrdVcor ):
            print "Convergence criteria met -- done."
            break

        #Log("PROGESS OF CORR.POTENTIAL FIT:\n")
       
        # how to do this general?
        vMfd = np.zeros((1))
        vMfd = MfdHam[0] - FSCoreH[0]
        #

        SkipDiis = dVsum > DiisThr and iMacroIt < DiisStart
        VcorLarge, dVcorLarge, vMfd, c0 = dc.Apply(VcorLarge, dVcorLarge, vMfd, Skip=SkipDiis)
        if not SkipDiis:
            print "\nVcor extrapolation: DIIS{{{:2d} {:2d}  {:8.2e}}}" %(dc.nDim, dc.iNext, c0)

        # how to do this general?
        MfdHam[0] = vMfd + FSCoreH[0]
        # 
        VcorLarge += dVcorLarge

    print "DMET SCF cycle converged after {} iterations" %(1+iMacroIt)
    print "Final residual on unit-cell vcor is {:.2e}"  %(dVsum)
