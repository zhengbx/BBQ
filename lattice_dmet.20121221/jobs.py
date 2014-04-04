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

"""support for scheduling lattice DMET calculations.
That is, tracking the dependencies of a calculation job on other jobs
(e.g., to provide initial guesses), and distributing them across MPI nodes."""

from copy import copy, deepcopy
from meanfield import FMeanFieldParams
from fragments import FDmetParams

def FmtParamsR(p, Prefix="p", Lines=None):
   if Lines is None:
      Lines = []
   Attr = p.__dict__
   Keys = Attr.keys()
   Keys.sort()
   for Key in Keys:
      Val = Attr[Key]
      if hasattr(Val,"__dict__"):
         FmtParamsR(Val,"%s.%s" % (Prefix, Key), Lines)
      else:
         Lines.append("{}.{} = {}".format(Prefix, Key, repr(Val)))
   return "\n".join(Lines)


class FJob(object):
   def __init__(self, Params, InputJobs = None):
      self.InputJobs = InputJobs  # if not None, take initial guess from here.
      self.InputJobs = [o for o in self.InputJobs if o is not None]
      if not self.InputJobs:
         self.InputJobs = None
      self.Params = Params
      pass

   def Run(self, Log, InputJobs):
      raise Exception("Run must be implemented!")

   def MarkAsDone(self):
      self.Done = True

   def __str__(self):
      return FmtParamsR(self.Params)

   def GetDependencies(self):
      #raise Exception("GetDependencies must be implemented!")
      if self.InputJobs is None:
         return None
      else:
         return self.InputJobs

class FDefaultParams(object):
   def __init__(self, LatticeWf):
      # lattice wave function: FWfDecl instance
      self.LatticeWf = LatticeWf
      # lattice mean-field settings
      self.MeanField = FMeanFieldParams()
      # dmet settings
      self.DMET = FDmetParams()
   def __str__(self):
      return FmtParamsR(self)

class FJobGroup(object):
   def __init__(self, Jobs):
      self.Jobs = Jobs

   def Run(self, Log, PrintJobSummary=None):
      """PrintJobSummary: if provided, a function called with a list
      of all jobs already done as argument."""

      # todo: do something more clever, like distrbuting stuff
      # over the network etc.
      JobsLeft = copy(self.Jobs)
      JobsDone = []

      while JobsLeft:
         DidSomething = False
         for Job in JobsLeft:
            Deps = Job.GetDependencies()
            #print "job:", Job
            #print "job.deps -> ", Deps
            #raise SystemExit
            if Deps is None or not [o for o in Deps if o not in JobsDone]:
               # prerequisites are met. Run current job.
               #try:
               if 1:
                  with Log.Section("j%02x"%(len(JobsDone)),"Invoking job '%s'" % Job.Params):
                     with Log.Section('par','Full parameter list'):
                        Log("%s"%Job)
                     Job.Run(Log, Deps)
                  Job.MarkAsDone()
               #except Exception,e:
               else:
                  Log(Log.WARNING, "job '{}' failed: {}", Job, e)
                  pass
               del JobsLeft[JobsLeft.index(Job)]
               JobsDone.append(Job)
               if PrintJobSummary is not None:
                  PrintJobSummary(Log, JobsDone)
               DidSomething = True
               break
               # ^- that's only here to get the job order more consistent with
               #    our input job order.
         if not DidSomething:
            Log(Log.WARNING, "FJobGroup: some jobs could not be scheduled, "
               "because they depend on jobs not controlled by this scheduler "
               "and prerequisites for running were never met.")
            break

      ## topologically sort job network by dependencies.
      #RootNodes = [o for o in self.Jobs if o.GetDependencies() is None]
      pass

def ToClass(kw):
   class _Struct:
      def __init__(self, kw):
         self.__dict__.update(kw)
   return _Struct(kw)
