# Copyright (c) 2012  Gerald Knizia
# 
# This file is part of the IR/WMME program
# (See http://www.theochem.uni-stuttgart.de/~knizia)
# 
# IR/WMME is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
# 
# IR/WMME is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/

from numpy import *
from scipy.linalg import *
import numpy as np
from diis import FDiisContext
from time import time

MakeIntegralsExecutable = "./wmme"
TmpDir = None
BasisLibDir = "./bases/"

set_printoptions(precision=8,linewidth=10060,suppress=True,threshold=nan)

pResultFmt = " %-32s%18.12f"
pResultFmtAnnoted = " %-32s%18.12f  (%s)"
pResultFmtS = " %-32s%18s"
pResultFmtI = " %-32s%5i"
pResultFmtIAnnoted = " %-32s%5i  (%s)"

THRORB = 1e-10      # threshold for orbital gradient in RHF
THRDEN = 1e-8       # threshold for change in energy in RHF
MAXIT_HF = 1002     # maximum num. iterations for RHF
ORB_OCC = 2.

ToAng =     0.5291772108  # molpro default.

def ElementNameDummy():
   ElementNames = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn".split()
   ElementNumbers = dict([(o,i+1) for (i,o) in enumerate(ElementNames)])
   return ElementNames, ElementNumbers
ElementNames, ElementNumbers = ElementNameDummy()

def mdot(*args):
   """chained matrix product: mdot(A,B,C,..) = A*B*C*...
   No attempt is made to optimize the contraction order."""
   r = args[0]
   for a in args[1:]:
      r = dot(r,a)
   return r

def dot2(A,B): return dot(A.flatten(),B.flatten())

class FAtomSet:
   def __init__(self, Positions, Elements, Orientations=None):
      """Positions: 3 x nAtom matrix. Given in atomic units (ABohr).
      Elements: element name (e.g., H) for each of the positions.
      Orientations: If given, a [3,3,N] array encoding the standard
      orientation of the given atoms (for replicating potentials!). For
      each atom there is a orthogonal 3x3 matrix denoting the ex,ey,ez
      directions."""
      self.Pos = Positions
      assert(self.Pos.shape[0] == 3 and self.Pos.shape[1] == len(Elements))
      self.Elements = Elements
      self.Orientations = Orientations
   def MakeXyz(self,NumFmt = "%15.8f"):
      Lines = []
      for i in range(len(self.Elements)):
         Lines.append(" %5s {0} {0} {0}".format(NumFmt) % (\
            self.Elements[i], self.Pos[0,i], self.Pos[1,i], self.Pos[2,i]))
      return "\n".join(Lines)
   def nElecNeutral(self):
      """return number of electrons present in the total system if neutral."""
      return sum([ElementNumbers[o] for o in self.Elements])
   def fCoreRepulsion(self):
      N = len(self.Elements)
      Charges = array([ElementNumbers[o] for o in self.Elements])
      fCoreEnergy = 0
      for i in xrange(N):
         for j in xrange(i):
            fCoreEnergy += Charges[i] * Charges[j] / norm(self.Pos[:,i] - self.Pos[:,j])
      return fCoreEnergy
   def __str__(self):         return self.MakeXyz()
   def __len__(self):         return len(self.Elements)
   def __getitem__(self,key): return self.Elements[key], self.Pos[:,key]

def MakeMolecularSystemRaw(Atoms, OrbBasis, FitBasis, BasisLibs):
   """setup an electronic system using the quantum chemistry Hamiltonian as
   defined by the given atoms and basis set names."""
   from os import path
   from tempfile import mkdtemp
   from shutil import rmtree
   from commands import getstatusoutput
   # make a directory to store our input/output in.
   BasePath = mkdtemp(prefix="bf-int.", dir=TmpDir)
   def Cleanup():
      rmtree(BasePath)

   # assemble arguments to integral generation program
   FileNameXyz = path.join(BasePath, "ATOMS")
   FileNameCoreH = path.join(BasePath, "INT1E")
   FileNameInt2e = path.join(BasePath, "INT2E")

   Args = []
   Args.append("--orb-trafo=Smh --matrix-format=npy")
   # ^- calculate integrals in symmetrically orthogonalized AO basis
   #    and store matrices in numpy .npy output matrix format.
   for BasisLib in BasisLibs:
      Args.append("--basis-lib=%s" % path.join(BasisLibDir, BasisLib))
   Args.append("--save-coreh=%s" % FileNameCoreH)
   Args.append("--save-fint2e=%s" % FileNameInt2e)
   Args.append("--basis-orb=%s" % OrbBasis)
   Args.append("--basis-fit=%s" % FitBasis)
   Args.append("--atoms-au=%s" % FileNameXyz)

   XyzLines = "%i\n\n%s\n" % (len(Atoms), Atoms.MakeXyz("%24.16f"))
   # ^- note on the .16f: it actually does make a difference. I had .8f
   #    there before, and it lead to energy changes on the order of 1e-8
   #    when treating only non-redundant subsystem out of a symmetric
   #    arrangement.
   try:
      with open(FileNameXyz, "w") as File:
         File.write(XyzLines)

      Cmd = "%s %s" % (MakeIntegralsExecutable, " ".join(Args))
      print "!Invoking %s\n" % Cmd
      iErr, Output = getstatusoutput(Cmd)
      if ( iErr != 0 ):
         raise Exception("Integral calculation failed. Output was:\n%s" % Output)

      CoreH = np.load(FileNameCoreH)
      Int2e = np.load(FileNameInt2e)
      print Output + "\n"
   except:
      Cleanup()
      raise
   # got everything we need. Delete the temporary directory.
   Cleanup()

   nOrb = CoreH.shape[0]
   Int2e = Int2e.reshape((Int2e.shape[0], nOrb, nOrb))
   CoreEnergy = Atoms.fCoreRepulsion()

   return CoreEnergy, CoreH, Int2e

def ReadXyzFile(FileName,Scale=1./ToAng):
   Text = open(FileName,"r").read()
   Lines = Text.splitlines()
   # allowed formats: <nAtoms> \n Desc \n <atom-list>
   #              or: <atom-list> (without any headers)
   # in the first case, only the first nAtoms+2 lines are read, in the
   # second case everything which does not look like a xyz line is
   # ignored.
   nAtoms = None
   r = 0,-1
   if ( len(Lines[0].split()) == 1 ):
      nAtoms = int(Lines[0].split()[0])
      r = 2,nAtoms+2
   Atoms = []
   Xyz = []
   for Line in Lines:
      ls = Line.split()
      try:
         Atom = ls[0]
         x,y,z = float(ls[1]), float(ls[2]), float(ls[3])
      except:
         continue
      Atom = Atom[0].upper() + Atom[1:].lower()
      # maybe we should allow for (and ignore) group numbers after the
      # elements?
      if Atom not in ElementNames:
         raise Exception("while reading '%s': unrecognized element '%s'." % (FileName,Atom))
      Atoms.append(Atom)
      Xyz.append((x,y,z))
   Xyz = Scale*array(Xyz).T
   if 0:
      print "*read '%s':\n%s" % (FileName, str(FAtomSet(Xyz, Atoms)))
   return Xyz, Atoms


def MakeFockOp(Rdm, CoreH, Int2e_Frs, Orb=None, ExchFactor=1.):
   """make closed-shell/spin ensemble Fock matrix:
         f[rs] = h[rs] + j[rs] - .5 k[rs].
   where j[rs] = (rs|tu) rdm[tu]
   and   k[rs] = (rs|tu) rdm[su]
   Return (f, j, k)."""
   if Orb is None:
      OccNum,Orb = eigh(-Rdm); OccNum *= -1.
   else:
      OccNum = diag(mdot(Orb.T, Rdm, Orb))

   # filter out occupation numbers slightly broken due to rounding.
   OccNum[OccNum < 0.] = 0.
   OccNum[OccNum > ORB_OCC] = ORB_OCC

   nFit, nOrb = Int2e_Frs.shape[0], Int2e_Frs.shape[1]

   # make coulomb.
   gamma = einsum("Frs,rs",Int2e_Frs, Rdm)
   j = einsum("Frs,F",Int2e_Frs,gamma)

   # make exchange -- first need to find the occupied orbitals.
   # (they can be non-standard, e.g., in atomic density initial guess)
   nOcc = len([o for o in OccNum if abs(o-ORB_OCC)<1e-8])
   OccOrb = Orb[:,:nOcc] * OccNum[:nOcc]**.5
   k1 = einsum("Frs,si->Fir", Int2e_Frs, OccOrb).reshape((nFit * nOcc, nOrb))
   k = dot(k1.T, k1)

   Fock = CoreH + j - (ExchFactor/ORB_OCC) * k
   return Fock, j, k


def RunHf(CoreEnergy, CoreH, Int2e_Frs, nOcc, Print):
   p = Print
   StartTime = time()
   nOrb = CoreH.shape[0]
   # coreH guess
   Fock = 1.*CoreH

   LastEnergy = 0.
   if p:
      print "*** HARTREE-FOCK\n"
      print " %-35s%5i FUNCTIONS" % ("ORBITAL BASIS:", nOrb)
      print " %-35s%5i FUNCTIONS" % ("FITTING BASIS:", Int2e_Frs.shape[0])
      print

   Converged = False
   dc = FDiisContext(12)
   if p: print "   {:^5s} {:^14s} {:^14s} {:^11s} {:>4s}  {:^11s}".format(
      "ITER.","ENERGY","ENERGY CHANGE", "GRAD", "DIIS","TIME")
   LastRdm = None
   for it in range(MAXIT_HF):
      # make new orbitals and new rdm.
      Ew,Orbs = eigh(Fock)

      Rdm = ORB_OCC * dot(Orbs[:,:nOcc], Orbs[:,:nOcc].T)

      # make closed-shell Fock matrix:
      #   f[rs] = h[rs] + j[rs] - .5 k[rs]
      # and HF electronic energy:
      #   E = dot(rdm,.5*(h+f))
      Fock, j, k = MakeFockOp(Rdm, CoreH, Int2e_Frs, Orb=Orbs)
      Energy = dot2(.5*Fock+.5*CoreH, Rdm) + CoreEnergy

      # calculate orbital gradient and energy
      OrbGrad = dot(Fock,Rdm)
      OrbGrad = OrbGrad - OrbGrad.T

      fOrbGrad = norm(OrbGrad.flatten())
      dEnergy = Energy - LastEnergy
      LastEnergy = Energy

      if p: print " {:5d} {:14.8f} {:+14.8f} {:11.2e} {:>6s} {:8.2f}".format(
         it+1, Energy, dEnergy, fOrbGrad, dc, time() - StartTime)

      if (fOrbGrad < THRORB and abs(dEnergy) < THRDEN):
         Converged = True
         break
      if (it == MAXIT_HF - 1):
         break # not converged.

      Fock, OrbGrad, c0 = dc.Apply(Fock, OrbGrad)

   if (not Converged and MAXIT_HF != 1):
      ErrMsg = "%s failed to converge."\
               "  NIT =%4i  GRAD=%8.2e  DEN=%8.2e" % (WF_TYPE, it+1, fOrbGrad, dEnergy)
      raise Exception(ErrMsg)
   elif p:
      print
      print pResultFmt % ("1e energy", dot2(Rdm, CoreH))
      print pResultFmt % ("2e energy", .5*dot2(Rdm, Fock - CoreH))
      if (CoreEnergy != 0.0):
         print
         print pResultFmt % ("Core energy", CoreEnergy)
         print pResultFmt % ("Electronic energy", Energy - CoreEnergy)
      print pResultFmt % ("Hartree-Fock energy", Energy)
      print

   return Energy

def TestMoleculesModule():
   FileName = "benzene.xyz"
   #FileName = "ch2.xyz"
   #FileName = "ne.xyz"
   Xyz,Elements = ReadXyzFile(FileName)
   Atoms = FAtomSet(Xyz, Elements)

   OrbBasis = "cc-pVDZ"
   #OrbBasis = "MINAO"
   FitBasis = "univ-JKFIT"
   BasisLibs = ["def2-nzvpp-jkfit.libmol", "minao.libmol", "emsl_cc-pVDZ.libmol"]
   CoreEnergy, CoreH, Int2e = MakeMolecularSystemRaw(Atoms, OrbBasis, FitBasis, BasisLibs)

   print "*** SYSTEM '%s' [%s/%s]\n" % (FileName, OrbBasis, FitBasis)

   RefEnergy = -230.722074147857
   HfEnergy = RunHf(CoreEnergy, CoreH, Int2e, nOcc=Atoms.nElecNeutral()/2, Print=True)
   if ( abs(HfEnergy - RefEnergy) > 1e-8 ):
      raise Exception("Something went wrong. HF energy doesn't match up with reference result."+\
                    "\nShould be:   %16.8f" % RefEnergy +\
                    "\nActually is: %16.8f (dE=%.8f)" % (HfEnergy,HfEnergy-RefEnergy))


if ( __name__ == "__main__" ):
   TestMoleculesModule()




