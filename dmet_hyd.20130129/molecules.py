# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Gerald Knizia and Garnet K.-L. Chan, 2012;
# (c) Princeton University, 2012

"""This module provides support for:
   - Representing the geometry and properties molecules/sets of atoms (FAtomSet)
   - Creating and modifying quantum systems (MakeMolecularSystem) representing
     molecular electronic systems (i.e., calculating their Hamiltonian and other
     properites based on the geometry and basis definition).
   - Keeping information about what the functions in a AO basis set represent
     (FBasisDesc)
"""
import itertools as it
itt = it
from os import path
from sys import argv
from copy import copy, deepcopy
from tempfile import mkdtemp
from numpy import *
import numpy as np
from numpy.random import random as uniform_random, seed
from scipy.linalg import *

from diis import FDiisContext
from helpers import *
import settings as g
from quantum_systems import *


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
   def fSubsetCoreRepulsion(self, AtomList):
      N = len(self.Elements)
      Charges = array([ElementNumbers[o] for o in self.Elements])
      fCoreEnergy = 0
      for i in AtomList:
         for j in xrange(N):
            if j == i: continue
            fCoreEnergy += .5 * Charges[i] * Charges[j] / norm(self.Pos[:,i] - self.Pos[:,j])
      return fCoreEnergy
   def __str__(self):
      return self.MakeXyz()
   def __len__(self):
      return len(self.Elements)
   def __iter__(self):
      for i in range(len(self)):
         yield self[i]
   def __getitem__(self,key):
      """return tuple AtomName, Position."""
      assert(type(key) is int)
      assert(key < len(self.Elements))
      return self.Elements[key], self.Pos[:,key]



def MakeMolecularSystemRaw(Atoms, WfDecl, OrbBasis, FitBasis, BasisLibs, BasisDesc):
   """setup an electronic system using the quantum chemistry Hamiltonian as
   defined by the given atoms and basis set names."""
   from os import path
   from shutil import rmtree
   from commands import getstatusoutput
   # make a directory to store our input/output in.
   BasePath = mkdtemp(prefix="int-pf.", dir=g.TmpDir)
   def Cleanup():
      rmtree(BasePath)

   CoreEnergy = Atoms.fCoreRepulsion()

   # assemble arguments to integral generation program
   FileNameXyz = path.join(BasePath, "ATOMS")
   FileNameCoreH = path.join(BasePath, "INT1E")
   FileNameInt2e = path.join(BasePath, "INT2E")

   Args = []
   Args.append("--orb-trafo=Smh --matrix-format=npy")
   # ^- calculate integrals in symmetrically orthogonalized AO basis
   #    and store matrices in numpy .npy output matrix format.
   for BasisLib in BasisLibs:
      Args.append("--basis-lib=%s" % path.join(g.BasisLibDir, BasisLib))
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
      WriteFile(FileNameXyz, XyzLines)

      Cmd = "%s %s" % (g.MakeIntegralsExecutable, " ".join(Args))
      #print "!%s" % Cmd
      iErr, Output = getstatusoutput(Cmd)
      if ( iErr != 0 ):
         raise Exception("Integral calculation failed. Output was:\n%s" % Output)

      CoreH = np.load(FileNameCoreH)
      Int2e = np.load(FileNameInt2e)
   except:
      Cleanup()
      raise
   # got everything we need. Delete the temporary directory.
   Cleanup()

   nOrb = CoreH.shape[0]
   Int2e = Int2e.reshape((Int2e.shape[0], nOrb, nOrb))
   if 0:
      for (i,j) in it.product(xrange(CoreH.shape[0]),xrange(CoreH.shape[1])):
         if ( abs(CoreH[i,j]) < 1e-10 ):
            CoreH[i,j] = 0
      for (i,j,k) in it.product(xrange(Int2e.shape[0]),xrange(Int2e.shape[1]),xrange(Int2e.shape[2])):
         if ( abs(Int2e[i,j,k]) < 1e-10 ):
            Int2e[i,j,k] = 0
   if 0:
      # explicitly symmetrize integrals. Should not be required.
      CoreH = 0.5 * (CoreH + CoreH.T)
      Int2e = .5 * (Int2e + Int2e.transpose([0,2,1]))

   CoreEnergy = Atoms.fCoreRepulsion()

   if g.WF_TYPE == "UHF":
      # split the integrals into spin-orbitals. even orbitals get
      # alpha spin, odd ones get beta spin.
      nOrb = CoreH.shape[0]
      nFit = Int2e.shape[0]
      CoreH = CombineSpinComps(CoreH, CoreH)
      Int2e_UHF = zeros((nFit,2*nOrb,2*nOrb))
      Int2e_UHF[:, ::2, ::2] = Int2e
      Int2e_UHF[:,1::2,1::2] = Int2e
      Int2e = Int2e_UHF

   return FSystem(CoreEnergy, CoreH, Int2e, WfDecl, BasisDesc=BasisDesc.ForAtoms(Atoms))

def AbsorbHfCore(FullSystem, Atoms, BasisDesc):
   """create a new molecular system with core orbitals and electrons removed
   (that is: completely. There will neither be basis functions nor electrons
    left for them)"""
   assert(BasisDesc is not None) # need to know which are the core orbitals.
   fs = FullSystem
   Hfr = fs.RunHf(Print = (g.IPRINT >= 1))
   nOrb = fs.nOrb

   # make a list of all core orbitals.
   FunctionDesc = BasisDesc.ForAtoms(Atoms)
   assert(len(FunctionDesc) == nOrb)
   CoreOrbitals = []
   for (i,FnDesc) in enumerate(FunctionDesc):
      if FnDesc.IsCore():
         CoreOrbitals.append(i)
   nCore = len(CoreOrbitals)
   if nCore == 0:
      return FullSystem
   if g.IPRINT >= 1:
      print g.pResultFmtS % ("Deleting core orbitals", " ".join(map(str,CoreOrbitals)))

   OtherOrbs = array([o for o in range(nOrb) if o not in CoreOrbitals])
   Orbs = Hfr.Orbs
   def FindCoreMos():
      if 0: return range(nCore)
      # ^- just take first nCore orbtials as core

      # project Hartree-Fock MOs onto (selected) core AOs in order to
      # figure out what MOs to delete. Typically those should be the
      # first nCore MOs (i.e., the nCore MOs with the lowest
      # energies), but if trying some unusual cores or in certain
      # ionic combinations, this might always be the case.
      CoreP_AO = 1. * zeros((nOrb,nOrb))
      for i in CoreOrbitals:
         CoreP_AO[i,i] = 1.
      CoreNess = diag(mdot(Orbs.T, CoreP_AO, Orbs))
      CoreNess = [(o,io) for (io,o) in enumerate(CoreNess)]
      CoreNess.sort(reverse=True)
      CoreMos = [io for (io,o) in enumerate(CoreNess[:nCore])]
      CoreMos.sort()
      maxvp, imaxvp = CoreNess[nCore]    # maximum core projection
      mincp, imincp = CoreNess[nCore-1]  # minimum valence projection
      if ( mincp/maxvp < 10. ):
         print "\nWARNING: Selected core AOs are not well separated at HF level." +\
               "\n         Minimum projection of core MOs onto core AOs:    %8.4f  (Ew = %12.6f)" % (mincp, Hfr.Ews[imincp]) +\
               "\n         Maximum projection of valence MOs onto core AOs: %8.4f  (Ew = %12.6f)" % (maxvp, Hfr.Ews[imaxvp]) +\
               "\n"
      if ( CoreMos != list(range(nCore)) ):
         print "\nWARNING: Core orbitals are not the first %i canonical HF orbitals!" +\
               "\n         Deleting MOs: %s" % "  ".join(["%6i" % io for (io,o) in CoreNess]) +\
               "\n         Core-ity:     %s" % "  ".join(["%6.3f" % o for (io,o) in CoreNess]) +\
               "\n"
      return array(CoreMos)
   CoreOrb = Orbs[:,FindCoreMos()]

   # projector onto the valence space.
   ValenceP = eye(nOrb) - dot(CoreOrb,CoreOrb.T)

   # take the vectors we expect to have non-zero overlap with the valence
   # basis. Note that we can orthogonalize at most N - nCore such vectors,
   # because nCore eigenvalues are zero. Thus, if we just take any such
   # vectors directly from (1-PCore)*id  (which is easy--we know what the
   # core is) and smh then we get exactly what we want.
   ValenceOrbs = ValenceP[:,OtherOrbs]
   Overlap = dot(ValenceOrbs.T, ValenceOrbs)
   ValenceOrbs = dot(ValenceOrbs, MakeSmh(Overlap))
   #print "(CoreOrb * CoreOrb.T) ValenceP2:\n", mdot(dot(CoreOrb,CoreOrb.T), ValenceOrbs)
   assert(norm(mdot(dot(CoreOrb,CoreOrb.T), ValenceOrbs).flatten()) < 1e-10)
   #print "Valence Orbitals:\n", ValenceOrbs

   CoreRdm = 2. * dot(CoreOrb, CoreOrb.T)
   CoreFock, j, k = MakeFockOp(CoreRdm, fs.CoreH, fs.Int2e_Frs, Orb=Orbs)
   DeltaCoreEnergy = .5*dot2(CoreFock + fs.CoreH, CoreRdm)
   print g.pResultFmt % ("Hartree-Fock energy", Hfr.Energy)
   print g.pResultFmt % ("Absorbed core energy", DeltaCoreEnergy)
   CoreEnergy = fs.CoreEnergy + DeltaCoreEnergy

   CoreH = mdot(ValenceOrbs.T, CoreFock, ValenceOrbs)
   #Int2e = einsum("Frs,rx,sy->Fxy", fs.Int2e_Frs, ValenceOrbs, ValenceOrbs)
   Int2e_1 = einsum("Frs,sy->Fry", fs.Int2e_Frs, ValenceOrbs)
   Int2e = einsum("Fry,rx->Fxy", Int2e_1, ValenceOrbs)
   del Int2e_1

   WfDecl = FWfDecl(fs.WfDecl.nElec - 2*nCore, fs.WfDecl.Ms2)
   BasisDesc1 = [fs.BasisDesc[o] for o in OtherOrbs]
   if 0:
      print "Reduced basis after core removal:"
      print " ".join(["   %4i   " % o for o in range(len(BasisDesc1))])
      print " ".join(["{:^10s}".format(o) for o in BasisDesc1])
   if g.IPRINT >= 1:
      print
   return FSystem(CoreEnergy, CoreH, Int2e, WfDecl, BasisDesc=BasisDesc1)

def MakeMolecularSystem(Atoms, WfDecl, OrbBasis, FitBasis, BasisLibs, BasisDesc=None):
   """Calculate integrals for full systems, then optionally run HF and delete
      the core orbitals and electrons from everything."""
   FullSystem = MakeMolecularSystemRaw(Atoms, WfDecl, OrbBasis, FitBasis, BasisLibs, BasisDesc)
   if not g.ABSORB_HF_CORE:
      return FullSystem
   else:
      return AbsorbHfCore(FullSystem, Atoms, BasisDesc)

def ReadXyzFile(FileName,Scale=1./ToAng):
   Text = ReadFile(FileName)
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

class FBasisFnDesc:
   def __init__(self, Type, Class, Element):
      from atom_configs import GetAtomAverageOccupation
      self.Type = Type
      self.Shell = int(Type[0])
      self.AngMom = Type[1]
      assert(self.AngMom in "spdfghi")
      assert(Class in ["Valence", "Core", "Pol"])
      self.Class = Class
      self.Element = Element
      self.iAtom = None
      self.AtomicOccupancy = GetAtomAverageOccupation(ElementNumbers[Element],
         self.Shell, "spdfghi".index(self.AngMom))
   def IsCore(self):
      return self.Class == "Core"
   def IsValence(self):
      return self.Class == "Valence"
   def iShell(self):
      return self.Shell
   def __repr__(self):
      Fmt = {"Valence": "%s",
             "Core": "[%s]",
             "Pol": "(%s)"}[self.Class]
      AtomId = ""
      if self.iAtom is not None:
         AtomId = "%2i " % self.iAtom
      return AtomId + self.Element + " " + Fmt % self.Type

class FBasisDesc(dict):
   def __init__(self,Functions=None):
      if Functions is not None:
         self.Functions = Functions
      pass
   def Add(self, Elements, Descs):
      """Elements: space-separated list of elements
         Desc: list of basis functions. Examples (F/cc-pVDZ):
               "[1s] 2s (3s) 2px 2py 2pz (3px) (3py) (3pz) (3d0) (3d2-) (3d1+) (3d2+) (3d1-)"
               "[1s] 2s (3s) 2px 2py 2pz (3px 3py 3pz) (3d0 3d2- 3d1+ 3d2+ 3d1-)"
               "[1s] 2s (3s) 2p (3p) (3d)
         (all three are equivalent)
      """
      Functions = []
      Status = "Valence"
      for s in Descs.split():
         if s.startswith("[") or s.startswith("("):
            assert(Status == "Valence")
            Status = {"[": "Core", "(": "Pol"}[s[0]]
            s = s[1:]
         NextStatus = Status
         if s.endswith("]") or s.endswith(")"):
            assert(Status == {"]": "Core", ")": "Pol"}[s[-1]])
            NextStatus = "Valence"
            s = s[:-1]
         Shell = int(s[0])
         Suffices = [""]
         if s.endswith("p"): Suffices = ["x","y","z"]
         if s.endswith("d"): Suffices = ["0","2-","1+","2+","1-"]
         for Suffix in Suffices:
            for Element in Elements.split():
               if Element not in self:
                  self[Element] = []
               self[Element].append( FBasisFnDesc(s + Suffix, Status, Element) )
         Status = NextStatus
         #print "%s <- %s" % (Element, Functions)
   def ForAtoms(self, AtomSet):
      Functions = []
      for (iAtom,Element) in enumerate(AtomSet.Elements):
         Functions += deepcopy(self[Element])
         for i in range(len(self[Element])):
            Functions[-i-1].iAtom = iAtom
      return Functions

def MakeDefaultBasisDesc(OrbBasis):
   BasisDesc = FBasisDesc()
   if OrbBasis in ["MINAO", "STO-6G"]:
      BasisDesc.Add("H", "1s")
      BasisDesc.Add("Li Be", "[1s] 2s")
      BasisDesc.Add("B C N O F Ne", "[1s] 2s 2px 2py 2pz")
      #print "\n\n!!!! ^- FIXME: remove 2s from C N F core.\n\n"
      BasisDesc.Add("Cl", "[1s 2s] 3s [2p] 3p")
      BasisDesc.Add("Br", "[1s 2s 3s] 4s [2p 3p] 4p [3d]")
      BasisDesc.Add("Fe", "[1s 2s 3s] (4s) (5s) [2p 3p] (4p) 3d")
   elif OrbBasis in ["cc-pVDZ"]:
      BasisDesc.Add("H", "1s (2s) (2p)")
      BasisDesc.Add("Li Be", "[1s] 2s (3s) (3px 3py 3pz)")
      BasisDesc.Add("B C N O F Ne", "[1s] 2s (3s) 2px 2py 2pz (3px 3py 3pz) (3d0 3d2- 3d1+ 3d2+ 3d1-)")
      #BasisDesc.Add("B C N O F Ne", "[1s] 2s (3s) 2p (3p 3d)")
   else:
      assert("basis not recognized -- need to know its composition to assign function types.")
   return BasisDesc


#def MakeIrrep0Basis(BasisDesc, Atoms):
   #def HasSymmetry(Op):
      #Pos1 = dot(Op,Atoms.Pos)
      #DistSq = ((Pos1.reshape(3,1,len(Atoms)) - Atoms.Pos.reshape(3,len(Atoms),1))**2).sum(axis=0)
      #for iAtom in range(len(Atoms)):
         #iOrig = argmin(DistSq[:,iAtom])
         #if ( DistSq[iOrig,iAtom] > 1e-10 or
              #Atoms.Elements[iOrig] != Atoms.Elements[iAtom] ):
            #return False
      #return True
      ##pass
   #def Mirror(xyz):
      #Op = eye(3)
      #if 'x' in xyz: Op[0,0] = -1
      #if 'y' in xyz: Op[1,1] = -1
      #if 'z' in xyz: Op[2,2] = -1
      #return Op
   #print HasSymmetry(Mirror('x'))
   #print HasSymmetry(Mirror('y'))
   #print HasSymmetry(Mirror('z'))
   #raise SystemExit
   #pass

def TestMoleculesModule():
   g.THRORB = 1e-10
   g.IPRINT = 2
   g.ABSORB_HF_CORE = True
   Xyz,Elements = ReadXyzFile("benzene.xyz")
   Atoms = FAtomSet(Xyz, Elements)

   WfDecl = FWfDecl(nElec=Atoms.nElecNeutral())
   OrbBasis = "cc-pVDZ"
   FitBasis = "univ-JKFIT"
   BasisLibs = [""]
   BasisDesc = MakeDefaultBasisDesc(OrbBasis)
   MakeIrrep0Basis(BasisDesc, Atoms)
   BasisLibs = ["def2-nzvpp-jkfit.libmol", "minao.libmol", "emsl_cc-pVDZ.libmol"]
   System = MakeMolecularSystem(Atoms, WfDecl, OrbBasis, FitBasis, BasisLibs, BasisDesc)
   RefEnergy = -230.722074147857
   hfr = System.RunHf(Print=True)
   if ( abs(hfr.Energy - RefEnergy) > 1e-8 ):
      raise Exception("Something went wrong. HF energy doesn't match up."+\
                    "\nShould be:   %16.8f" % RefEnergy +\
                    "\nActually is: %16.8f (dE=%.8f)" % (hfr.Energy,hfr.Energy-RefEnergy))


if ( __name__ == "__main__" ):
   TestMoleculesModule()




