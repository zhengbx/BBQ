# This file is part of the lattice-dmet program. lattice-dmet is free
# software: you can redistribute it and/or modify it under the terms of
# the GNU General Public License as published by the Free Software
# Foundation, version 3.
# 
# lattice-dmet is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with bfint (LICENSE). If not, see http://www.gnu.org/licenses/
# 
# Authors:
#    Gerald Knizia, 2012

from output import *
from textwrap import dedent
from meanfield import *
from copy import copy
from lattice_model import FLatticeModel
from fragments import FFragment, FDmetContext, FDmetParams
from jobs import FJob, FJobGroup, FDefaultParams, ToClass

Banner = """\
___________________________________________________________________

    L A T T I C E   D M E T                             [v20121221]
_ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
"""


class FHubbardModel2dHoneyComb(FLatticeModel):
   def __init__(self, t, U):
      UnitCell = []
      UnitCell.append( ("U1", np.array([0.,0.])) )
      UnitCell.append( ("U2", np.array([1.,0.])) )
      self.t = t
      self.U = U

      ToRad = 2.*np.pi/180.
      c = -np.cos(120. * ToRad)
      s =  np.sin(120. * ToRad)
      #print c, s
      #raise SystemExit

      #
      #    1-----2
      #           \
      #            3-----4
      #           /
      #    5-----6
      #
      # we use as lattice translations the vectors:
      #    1 -> 5 (that's ``down by 2*sin(120 deg)'' and
      #    1 -> 3 (that's ``right by 1-cos(120 deg) and down by sin(120 deg)'')

      LatticeVectors = np.zeros((2,2),float)
      LatticeVectors[0,:] = [  0., 2*s]
      LatticeVectors[1,:] = [1.+c,   s]
      LatticeVectors = LatticeVectors.T
      #print LatticeVectors
      #raise SystemExit

      # FIXME: MaxRangeT handling is wrong for tilted lattices
      # in super-cell code.
      FLatticeModel.__init__(self, UnitCell, LatticeVectors,
         ["U", "t"], MaxRangeT=1.5, EnergyFactor=1./len(UnitCell))
   def GetTij(self, (SiteTypeI,XyzI), (SiteTypeJ,XyzJ)):
      dXyz = XyzJ - XyzI
      if sum(dXyz**2) <= 1.01: return -self.t
      return 0.
   def MakeTijMatrix(self, SitesR, SitesC):
      # This function is technically unnecessay, it just calculates the tij
      # matrix in matrix form instead of element-wise as the default
      # implementation. Its presence will, however, significantly speed up
      # calculations on large lattices.
      tij = np.zeros((len(SitesR), len(SitesC)))
      dXyz = SitesR.Xyzs[:,np.newaxis,:] - SitesC.Xyzs[np.newaxis,:,:]
      tij[np.sum(dXyz**2,axis=2) <= 1.01] = -self.t
      return tij
   def GetUi(self, (SiteTypeI,XyzI)):
      return self.U

class FParams1(FDefaultParams):
   def __init__(self, Model, SuperCell, LatticeWf, Fragments):
      #FDefaultParams.__init__(self, Model, SuperCell, LatticeWf, Fragments)
      FDefaultParams.__init__(self, LatticeWf)
      self.Model = Model
      self.SuperCell = SuperCell
      self.Fragments = Fragments

      self.MeanField.MaxIt = 40
      self.DMET.UseInt2e = False
      #self.DMET.DiisThr = 1e99
      #self.DMET.DiisStart = 0
   def __str__(self):
      return "U = %10.5f  <n> = %10.5f" % (self.Model.U,
      self.LatticeWf.nElec/(1.*np.prod(self.SuperCell.TotalSize)))

class FHub1dJob(FJob):
   def __init__(self, Params, InputJob = None):
      # Input job: nothing or job to take initial guess for
      # mean field/vloc/embedded wave functions from.
      assert(InputJob is None or isinstance(InputJob, FHub1dJob))
      FJob.__init__(self, Params, [InputJob])
      pass
   def Run(self, Log, InputJobs):
      P = self.Params

      StartingGuess = {}
      if InputJobs is not None:
         assert(len(InputJobs) == 1)
         GuessJob = InputJobs[0]
         StartingGuess = GuessJob.Results
      else:
         StartingGuess['fock'] = P.InitialGuessSpinBias
         if P.MeanField.MaxIt == 1 and P.InitialGuessSpinBias is not None:
            StartingGuess['vcor'] = np.diag(P.InitialGuessSpinBias)

      sc = FSuperCell(P.Model, OrbType=P.LatticeWf.OrbType, **P.SuperCell.__dict__)
      LatticeSystem = FLatticeSystem(WfDecl=P.LatticeWf, SuperCell=sc,
         InitialGuess=StartingGuess['fock'], Params=P.MeanField, Log=Log)

      LatticeSystem.RunHf(Log=Log)
      DmetContext = FDmetContext(LatticeSystem, P.Fragments, P.DMET, Log)
      DmetResult = DmetContext.Run(Log, StartingGuess)

      # next step: how to organize DMET and the fragmentation,
      # and how to organize storage of results.
      # (need something for vcor, fock, and potentially CI vectors)
      # this.. is probably not it. Just to get stuff running.
      # Note in particular that this just keeps on accumulating Fock
      # matrices in memory which never get released until all job
      # objects are destroyed.
      self.Results = {
         'vcor': DmetResult.FullSystemVcor,
         'dmet': DmetResult,
         'fock': DmetResult.FullSystemFock,
         'Mu':  LatticeSystem.Mu,
         'Gap': LatticeSystem.Gap,
         'nSites': len(sc.Sites), # number of sites in the super-cell
      }
   def GetResultTable(self):
      P = self.Params
      R = self.Results
      Out = {}
      if P.LatticeWf.OrbType == "UHF":
         # in the UHF case, we get band gaps and chemical potentials
         # separately for alpha and beta spin. These are supplied as
         # tuples in LatticeSystem.Gap/.Mu. Unpack them.
         Out = {
            "Gap[A]": R["Gap"][0],
            "Gap[B]": R["Gap"][1],
            "Mu[A]": R["Mu"][0],
            "Mu[B]": R["Mu"][1],
            '<Sz>': (P.LatticeWf.nElecA() - P.LatticeWf.nElecB()) / (1.*R["nSites"]),
         }
      else:
         assert(P.LatticeWf.OrbType == "RHF")
         Out = {
            "Mu": R["Mu"],
            "Gap": R["Gap"],
         }
      Out.update({
         "U": P.Model.U,
         "Fragments": P.Fragments,
         "E/Site": R["dmet"].TotalEnergy,
         "<n>[DMET]": R["dmet"].TotalElec,
         '<n>': P.LatticeWf.nElec / (1.*R["nSites"]),
         "nSites": R["nSites"]*int(P.LatticeWf.OrbOcc())/2
      })
      return Out, {
         "<Sz>": "average number of spin-up minus spin-down electrons per site",
         "<n>": "charge density per site at mean-field level (as input)",
         "<n>/DMET": "charge density per site at DMET level",
      }

def PrintTable(Log, Results, Desc = None, SortKeys = None):
   if SortKeys is None: SortKeys = []

   def CmpResults(a,b):
      for Key in SortKeys:
         ic = cmp(Key in a, Key in b)
         if ic != 0: return ic
         if Key not in a: continue # neither has this result. go on.
         ic = cmp(a[Key], b[Key])
         if ic != 0:
            return ic
      return 0
   Results.sort(CmpResults)

   # make a list of all keys, for the table captions.
   AllKeys = set([])
   for Result in Results:
      AllKeys |= set(Result.keys())
   def CmpSortKeys(a,b):
      ic = cmp(a in SortKeys, b in SortKeys)
      if ic != 0: return -ic # keys which we sort by go first.
      if a in SortKeys:
         # primary sorting keys go first.
         return cmp(SortKeys.index(a), SortKeys.index(b))
      else:
         # sort rest lexicographically, by name
         return cmp(a,b)
   AllKeys = list(AllKeys)
   AllKeys.sort(CmpSortKeys)

   def FmtV(v):
      if isinstance(v,float):
         return "%14.8f" % v
      return str(v).replace(" ","")
      # ^- replace: to simplify parsing the result tables. space only
      #             used as separator between columns.

   Lines = []
   for Result in Results:
      Line = []
      for Key in AllKeys:
         if Key not in Result:
            Line.append("")
         else:
            v = Result[Key]
            Line.append(FmtV(v))
      Lines.append(Line)

   Fmts = []
   for iCol in range(len(AllKeys)):
      Max = len(AllKeys[iCol])
      for Line in Lines:
         Max = max(Max, len(Line[iCol]))
      Fmts.append("{:^%is}" % Max)

   Caption = "  ".join(Fmt.format(Key) for (Key,Fmt) in zip(AllKeys,Fmts))
   Log(Caption + "\n" + "-" * len(Caption))
   for Line in Lines:
      Log("  ".join(Fmt.format(Val) for (Val,Fmt) in zip(Line,Fmts)))


   if Desc:
      Log()
      Log("Notes:")
      for (ItemName, ItemDesc) in Desc.items():
         Log("  {:18} {}", ItemName, ItemDesc)

def main():
   np.set_printoptions(precision=3,linewidth=10060,suppress=False,threshold=np.nan)
   Log = FOutputLog()
   Log(Banner)

   Jobs = []
   for nImp in [4]:
      Fragments = [('FCI',list(range(nImp)))]
      #Fragments = [('CI(2)',list(range(nImp)))]
      if 1:
         FModelClass = FHubbardModel2dTilted
         FModelClass = FHubbardModel2dHoneyComb
         assert(nImp==4) # needs to be translated into ClusterSize. Note that
                         # cluster size comes in terms of unit cells, each of
                         # which has two sites in these models.
         ScParams = ToClass({'TotalSize': [16,16], 'PhaseShift': [-1,-1], 'ClusterSize': [2,1]})
      else:
         FModelClass = FHubbardModel1d
         ScParams = ToClass({'TotalSize': [480], 'PhaseShift': [-1], 'ClusterSize': [nImp]})
      ModelParams = {'t': 1., "U": 1.0 }
      Model = FModelClass(**ModelParams)

      # make a super-cell in order to determine fillings with unique ground
      # states at tight-binding level.
      SuperCell = FSuperCell(Model, **ScParams.__dict__)
      AllowedOccupations = SuperCell.CalcNonDegenerateOccupations(ThrDeg=1e-5)
      Occs = AllowedOccupations[:len(AllowedOccupations)/2+1]
      StartGuessNextU = None

      Occs = Occs[-1:] # <- only half filling
      for U in [1.,2.,3.,4.,5.,6.]:
      #for U in [1.]:
         StartGuess = StartGuessNextU # <- start with last U's half-filling result.
         StartGuessNextU = None
         ModelParams = {'t': 1., "U": U }
         Model = FModelClass(**ModelParams)

         for Occ in reversed(Occs):
            sp = 0
            #sp = 1 # spin polarization.
            iOccA = AllowedOccupations.index(Occ)
            LatticeWf = FWfDecl(nElecA=AllowedOccupations[iOccA+sp],
                              nElecB=AllowedOccupations[iOccA-sp],
                              OrbType="UHF")
            P = FParams1(
               Model=Model,#ToClass(ModelParams),
               SuperCell=ScParams,
               LatticeWf=LatticeWf,
               Fragments=Fragments)
            #P.DMET.DiisStart = 0
            #P.DMET.DiisThr = 1e-40
            P.DMET.MaxIt = 40
            P.MeanField.MaxIt = 1 # disable iterations.
            P.MeanField.DiisStart = 8
            P.InitialGuessSpinBias = None
            if 1:
               # add some bias to make the system preferrably go into
               # anti-ferromagnetic solutions. Otherwise we will always
               # get RHF solutions via UHF, even if the UHF solution is
               # lower.
               bias = ModelParams["U"]/2
               shift = 0.
               shift = ModelParams["U"]/2 if P.MeanField.MaxIt == 1 else 0.
               if P.LatticeWf.OrbType == "UHF":
                  P.InitialGuessSpinBias = np.zeros(2*len(SuperCell.UnitCell))
                  P.InitialGuessSpinBias[ ::4] += shift + bias
                  P.InitialGuessSpinBias[1::4] += shift - bias
                  P.InitialGuessSpinBias[2::4] += shift - bias
                  P.InitialGuessSpinBias[3::4] += shift + bias
                  #print P.InitialGuessSpinBias
                  #raise SystemExit
               else:
                  P.InitialGuessSpinBias = np.zeros(1*len(SuperCell.UnitCell))
                  P.InitialGuessSpinBias[:] = shift

            Jobs.append(FHub1dJob(P, InputJob=StartGuess))
            StartGuess = Jobs[-1]
            if 0:
               # enable this to propagate starting guesses from one U
               # to another. Here disabled for the honeycomb hubbard case,
               # because at some U between 3. and 4. the solution changes character.
               if StartGuessNextU is None:
                  StartGuessNextU = Jobs[-1]
            #if len(Jobs) == 1: break

   def PrintJobResults(Log, JobsDone):
      with Log.Section("r%03x" % len(JobsDone), "RESULTS AFTER %i OF %i JOBS:" % (len(JobsDone), len(Jobs)), 1):
         AllResults = []
         for Job in JobsDone:
            Result, Desc = Job.GetResultTable()
            AllResults.append(Result)
         PrintTable(Log, AllResults, Desc, SortKeys=["U", "<n>", "Fragments"])

   JobGroup = FJobGroup(Jobs)
   JobGroup.Run(Log, PrintJobResults)


main()

