# Draft for the main driver routine for joint dmet code.
# It needs the classes: HAM, GEOM, TYPE scetched out below.
# Here, TYPE is variable that is assigned any of the allowed
# TYPE-classes: Normal, BCS, Eph.
# No print-out has been inserted yet.
# It is not yet connected with the input module.
import sys

import numpy as np
import itertools as it

from helpers import mdot, ExtractSpinComp, CombineSpinComps
from vcor_fit import FitVcorComponent
from results import FDmetResult
from diis import FDiisContext
from inputs import Input
from normalDmet import NormalDmet
from geometry import BuildLatticeFromInput, Wavefct
from ham import Hamiltonian


def FitCorrelationPotential(Input, GEOM, TYPE, EmbMfdHam, nEmb, RdmHl):
    """fit a matrix nEmb x nEmb, referring to an operator connecting  
    all embedded sites to each other, to more closely approach the 
    correlated calculation result in with the mean-field method.
    Returns nEmb x nEmb matrix (always), even if due to the fitting
    criterion used the actual matrix is only a subset of this."""

    nImp = GEOM.nImp
    OrbType = GEOM.OrbType
    ORB_OCC = GEOM.OrbOcc()
    nElec = GEOM.nElec
    VcorType = Input.FITTING.VcorType
    VcorFitType = Input.FITTING.VcorFitType

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



def main(inputdict):

    Inp = Input(inputdict)

    Lattice = BuildLatticeFromInput(Inp.GEOMETRY, Inp.WAVEFUNCTION.OrbType)
    HAM = Hamiltonian(Inp.HAMILTONIAN)
    WAVEFCT = Wavefct(Inp.WAVEFUNCTION, Lattice)
    TYPE = NormalDmet(Inp.MFD, WAVEFCT.OrbType, WAVEFCT.nElec, WAVEFCT.nElecA, WAVEFCT.nElecB, WAVEFCT.Ms)
  
    #print Lattice.UnitCell print function not implemented yet
    Lattice.set_Hamiltonian(HAM)

    DmetMaxIt = Inp.DMET.max_iter
    ThrdVcor = Inp.DMET.conv_threshold
    DiisThr = Inp.DMET.diis_thr
    DiisStart = Inp.DMET.diis_start
    DiisDim = Inp.DMET.diis_dim


    MfdHam = Lattice.get_h0()
    FSCoreHam = Lattice.get_h0()
    Int2e = HAM.get_Int2e(Lattice, WAVEFCT.OrbType)

    Fragments = Lattice.supercell.fragments
  
    dc = FDiisContext(DiisDim)
    VcorLarge = TYPE.GuessVcor(WAVEFCT.OrbType, Inp.DMET, Lattice.supercell, HAM.U)
 
    for iMacroIt in range(DmetMaxIt):
        MfdHam_aug = np.zeros(MfdHam.shape)
        MfdHam_aug += MfdHam
        MfdHam_aug[0] += VcorLarge
        if HAM.transinv == True:
            # use translational symmetry
            MfdHamK_aug = Lattice.FFTtoK(MfdHam_aug) 
            MfdResult = TYPE.RunMfd(MfdHamK_aug, Int2e, Lattice, HamBlockDiag = True) 
            MfdRdm = Lattice.FFTtoT(MfdResult.Rdm)
        else:
            MfdResult = TYPE.RunMfd(MfdHam_aug, HamBlockDiag = False) 
            MfdRdm = MfdResult.Rdm
        FragmentResults = []
        FragmentPotentials = []
        for (iFragment,Fragment) in enumerate(Fragments):
            if Fragment.get_emb_method() is not None:
                EmbBasis = TYPE.MakeEmbBasis(MfdRdm.reshape(
                                             (MfdRdm.shape[0]*MfdRdm.shape[1],) 
                                             + MfdRdm.shape[2:]), 
                                             Fragment.get_sites())
                # check EmbBasis correct
                # trafo Int2e, n2eOrb
                EmbHam, EmbMfdHam, Int2e, EmbRdm = TYPE.MakeEmbHam(Lattice, 
                                                                  MfdHam_aug, MfdRdm, 
                                                                  HAM, EmbBasis, 
                                                                  Fragment.get_sites())
                HlResults = TYPE.ImpSolver(EmbCoreH, EmbMfdHam, EmbRdm, 
                                          Fragment.get_emb_method(), 
                                          Fragment.get_sites(), Int2e)
                FragmentResults.append((Fragment, HlResults))
                vloc = FitCorrelationPotential(Inp, GEOM, TYPE, EmbMfdHam, 
                                              TYPE.MakeEmbBasis.nEmb, RdmHl)
                FragmentPotentials.append(vloc)
        ClusterFactor = GEOM._EnergyFactor()
        # ^- for normalization with super-cell cluster size. Purely cosmetic.
        FragRes = FDmetResults(FragmentResults, ClusterFactor)
        FragRes.MeanFieldEnergy = MfdResult.MeanFieldEnergy
        FragRes.TotalEnergy = FragRes.TotalElecEnergy + MfdResult.CoreEnergy
        FragRes.FragmentPotentials = FragmentPotentials
        dVcorLarge = np.zeros(VcorLarge.shape)

        for (Fragment, VcorFull) in zip(GEOM.Fragments, FragRes.FragmentPotentials):
            nImp = len(Fragment.Sites)
            Vcor = Fragment.Factor * VcorFull[:nImp,:nImp]
            for Replica in Fragment.Replicas:
                VcorR = mdot(Replica.Trafo, Vcor, Replica.Trafo.T)
                E = list(enumerate(Replica.Sites))
                for (i,iSite),(j,jSite) in it.product(E,E):
                    dVcorLarge[iSite,jSite] += VcorR[i,j]

        dVsum = np.dot(dVcorLarge.flatten(),dVcorLarge.flatten())

        if (iMacroIt == 0 ):
            InitialHfEnergy = FSMfdResult.MeanFieldEnergy
        FragRes.dVc = dVsum
        if ( dVsum < ThrdVcor ):
            print "Convergence criteria met -- done."
            break

        #Log("PROGESS OF CORR.POTENTIAL FIT:\n")

        vMfd = np.zeros((1))
        vMfd = MfdHam[0] - FSCoreHam[0]
        SkipDiis = dVsum > DiisThr and iMacroIt < DiisStart
        VcorLarge, dVcorLarge, vMfd, c0 = dc.Apply(VcorLarge, dVcorLarge, vMfd, Skip=SkipDiis)
        if not SkipDiis:
            print "\nVcor extrapolation: DIIS{{{:2d} {:2d}  {:8.2e}}}" %(dc.nDim, dc.iNext, c0)
        MfdHam[0] = vMfd + FSCoreHam[0]
        VcorLarge += dVcorLarge
        
    print "DMET SCF cycle converged after {} iterations" %(1+iMacroIt)
    print "Final residual on unit-cell vcor is {:.2e}"  %(dVsum)


def parseOpt():
    if(len(sys.argv) != 2):
        raise Exception("No input file or more than one input file specified.")
    else:
        np.set_printoptions(precision=3,linewidth=10060,suppress=False,threshold=np.nan)
    inputfile = sys.argv[1]
    moreoptions = sys.argv[1:]
    return inputfile, moreoptions

if __name__ == '__main__':
    #inpfile = parseOpt()[0]
    sites = [(np.array([0., 0.]), "X")]
    shape = np.array([
        [1., 0.],
        [0., 1.],
    ])
    inpdic = {
        'HAMILTONIAN': {'Type': 'Hubbard', 'U': 3},
        'WAVEFUNCTION': {'OrbType': 'UHF', 'filling': 0.5, 'Ms': 0},
        'GEOMETRY':
        {'UnitCell': {'Sites':sites, 'Shape':shape},
         'ClusterSize': np.array([2, 2]),
         'LatticeSize': np.array([8, 8]),
         'Fragments': [{"Sites":range(4), "ImpSolver":'Fci',
                        "Fitting":'FullRdm'},],
         'BoundaryCondition': 'pbc',},
       'DMET':
       {'init_guess_type': 'AF'},
       'MFD':
       {'mfd_algo': 'hf'}
    }
    main(inpdic)

