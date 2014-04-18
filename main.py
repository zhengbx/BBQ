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
from results import FDmetResult
from diis import FDiisContext
from inputs import Input
from normalDmet import NormalDmet
from geometry import BuildLatticeFromInput, Wavefct
from ham import Hamiltonian





def main(inputdict):

    Inp = Input(inputdict)
    print Inp    

    Lattice = BuildLatticeFromInput(Inp.GEOMETRY, Inp.WAVEFUNCTION.OrbType)
    HAM = Hamiltonian(Inp.HAMILTONIAN)
    WAVEFCT = Wavefct(Inp.WAVEFUNCTION, Lattice)
    TYPE = NormalDmet(Inp.MFD, Inp.FITTING, WAVEFCT.OrbType, WAVEFCT.nElec, WAVEFCT.nElecA, WAVEFCT.nElecB, WAVEFCT.Ms)
  
    #print Lattice.UnitCell print function not implemented yet
    Lattice.set_Hamiltonian(HAM)

    DmetMaxIt = Inp.DMET.max_iter
    ThrdVcor = Inp.DMET.conv_threshold
    DiisThr = Inp.DMET.diis_thr
    DiisStart = Inp.DMET.diis_start
    DiisDim = Inp.DMET.diis_dim


    MfdHam = Lattice.get_h0()
    CoreH = Lattice.get_h0()
    Int2e = HAM.get_Int2e(Lattice, WAVEFCT.OrbType)

    Fragments = Lattice.supercell.fragments
  
    dc = FDiisContext(DiisDim)
    VcorLarge = TYPE.GuessVcor(WAVEFCT.OrbType, Inp.DMET, Lattice.supercell, HAM.U)
 
    for iMacroIt in range(DmetMaxIt):

        MfdHam_aug = np.zeros(MfdHam.shape)
        MfdHam_aug += MfdHam
        MfdHam_aug[0] += VcorLarge

        CoreH_aug = np.zeros(CoreH.shape)
        CoreH_aug += CoreH
        CoreH_aug[0] += VcorLarge
               
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
                EmbCoreH, EmbCoreH_aug, EmbInt2e, EmbRdm = TYPE.MakeEmbHam(Lattice, 
                                                                          MfdHam_aug, MfdHam, 
                                                                          MfdRdm, Int2e, EmbBasis, 
                                                                          Fragment.get_sites())
                HlResults = TYPE.ImpSolver(EmbCoreH, EmbRdm, 
                                          Fragment.get_emb_method(), 
                                          Fragment.get_sites(), EmbInt2e)
                FragmentResults.append((Fragment, HlResults))
                vloc = TYPE.FitCorrelationPotential(Inp.FITTING, Fragment.get_sites(), 
                                                   HlResults.nEmbElec, EmbCoreH_aug, 
                                                   HlResults.Rdm)
                FragmentPotentials.append(vloc)
        FragRes = FDmetResult(FragmentResults)#, ClusterFactor)
        FragRes.MfdEnergy = MfdResult.Energy
        FragRes.TotalEnergy = FragRes.TotalElecEnergy + HAM.get_core_energy(Lattice)
        FragRes.FragmentPotentials = FragmentPotentials

        dVcorLarge = TYPE.VcorUpdate(FragmentPotentials, Fragments, VcorLarge.shape)
        dVsum = np.dot(dVcorLarge.flatten(),dVcorLarge.flatten())

        if (iMacroIt == 0 ):
            InitialHfEnergy = MfdResult.Energy
        FragRes.dVc = dVsum
        if ( dVsum < ThrdVcor ):
            print "Convergence criteria met -- done."
            break

        #Log("PROGESS OF CORR.POTENTIAL FIT:\n")
        vMfd = np.zeros((1))
        vMfd = MfdHam[0] - CoreH[0]
        SkipDiis = dVsum > DiisThr and iMacroIt < DiisStart
        VcorLarge, dVcorLarge, vMfd, c0 = dc.Apply(VcorLarge, dVcorLarge, vMfd, Skip=SkipDiis)
        if not SkipDiis:
            print "\nVcor extrapolation: DIIS{{{:2d} {:2d}  {:8.2e}}}" %(dc.nDim, dc.iNext, c0)
        MfdHam_aug[0] = vMfd + CoreH_aug[0]
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
                        "Fitting":'FullRdm', "Factor": 1},],
         'BoundaryCondition': 'pbc',},
       'DMET':
       {'init_guess_type': 'AF'},
       'MFD':
       {'mfd_algo': 'diag1el'},
       'FITTING':
       {'v_fit_domain': 2, 'vfit_method': 'FullRdm'}
    }
    main(inpdic)

