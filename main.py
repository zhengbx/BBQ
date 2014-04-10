# Draft for the main driver routine for joint dmet code.
# It needs the classes: HAM, GEOM, TYPE scetched out below.
# Here, TYPE is variable that is assigned any of the allowed
# TYPE-classes: Normal, BCS, Eph.
# No print-out has been inserted yet.
# It is not yet connected with the input module.
import sys

import numpy as np
import itertools as it
import ast

from helpers import mdot, ExtractSpinComp, CombineSpinComps
from vcor_fit import FitVcorComponent
from results import FDmetResult
from diis import FDiisContext
from input import Input
from normalDmet import NormalDmet
from geometry import BuildLatticeFromInput
# .. work in progress
#from hamiltonian import Hamiltonian, Geometry
#from lattice_model import *

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



def main(argv):

   if(len(argv) !=1):
      raise Exception("No input file or more than one input file specified.")
   else:
      np.set_printoptions(precision=3,linewidth=10060,suppress=False,threshold=np.nan)

   #input parameters should be read from file
   inputfile = open(argv[0], 'r')
   obj_code = 'dict({})'.format(inputfile.read())
   inputdict = eval(obj_code)
   Inp = Input(inputdict)
   #Inp = Input({'DMET':{'max_iter':10},'MFD':{'scf_solver':'UHF'}})

   Lattice = BuildLatticeFromInput(Inp.GEOMETRY)
   
   #print Lattice.UnitCell print function not implemented yet
   HAM = Hamiltonian(Inp.HAMILTONIAN)

   Lattice.set_Hamiltonian(HAM)

   DmetMaxIt = Inp.DMET.max_iter
   ThrdVcor = Inp.DMET.conv_threshold
   DiisThr = Inp.DMET.diis_thr
   DiisStart = Inp.DMET.diis_start
   DiisDim = Inp.DMET.diis_dim

   TYPE = NormalDmet

   MfdHam = Lattice.get_h0()
   MfdHam += initguess

   FSCoreHam = Lattice.get_h0()
  
   Fragments = Lattice.supercell.fragments

   dc = FDiisContext(DiisDim)

   if 'vcor' in StartingGuess:
       VcorLarge = 1.*StartingGuess['vcor']
   else:
       #replace with data from GEOM
       VcorLarge = np.zeros_like(MfdHam[0,:,:])


   for iMacroIt in range(DmetMaxIt):
       MfdHam_aug = MfdHam + VcorLarge
       MfdResult = TYPE.RunMfd(MfdHam_aug)
       Rdm = MfdResult.Rdm
       FragmentResults = []
       FragmentPotentials = []
       for (iFragment,Fragment) in enumerate(Fragments):
         if Fragment.get_emb_method() is not None:
           EmbBasis = TYPE.MakeEmbBasis(Rdm, Fragment.get_sites())
           EmbHam, EmbMfdHam = TYPE.MakeEmbHam(EmbBasis, MfdHam, HAM, Fragment.get_sites())
           HlResults = TYPE.ImpSolver(EmbHam, EmbMfdHam, Fragment.get_emb_method())
           FragmentResults.append((Fragment, HlResults))
           vloc = FitCorrelationPotential(Inp, GEOM, TYPE, EmbMfdHam, TYPE.MakeEmbBasis.nEmb, RdmHl)
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


main(sys.argv[1:])

