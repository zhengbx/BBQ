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

"""Support for hierarchical recording the progress of a calculation group, and
for finding out which calculations were done where, and with which results.
(note: the python logging module is nice and all, but does not do what we want!)"""
from sys import stdout, stderr
import numpy as np

pResultFmtF = "%-34s%14.8f"
pResultFmtS = "%-34s%14s"
pResultFmtI = "%-34s%14i"

def FmtMat(m,ind=""):
   Lines = []
   def FmtFlt(o):
      w,p = 10,5
      BaseFmt = "%{}.{}f".format(w,p)
      s = BaseFmt % o
      while s[-1] == '0':
         s = s[:-1]
      s += (w-len(s))*" "
      return s
   for Line in m:
      Lines.append(ind + " [" + " ".join(FmtFlt(o) for o in Line) + "]")
   Lines[0] = Lines[0].replace(" [","[[")
   Lines[-1] += "]"
   return "\n".join(Lines)

def pResultFmt(Name, obj, Annotation=None):
   if isinstance(obj, (float,np.float64)):
      r = pResultFmtF % (Name,obj)
   elif isinstance(obj, (complex,np.complex128)):
      r = "{0:<34s}{1.real:14.8f}{1.imag:+14.8f}j".format(Name,obj)
   elif isinstance(obj, (int,np.int32,np.int64)):
      r = pResultFmtI % (Name,obj)
   elif isinstance(obj, np.ndarray):
      if obj.shape > 1:
         if Annotation is None:
            r = "%-34s\n%s" % (Name+':',FmtMat(obj,ind="    "))
         elif Annotation is "UHF":
            a,b = obj[::2,::2], obj[1::2,1::2]
            r =    "%-34s\n%s" % (Name+' (charge):',FmtMat(a+b,ind="    "))
            r += "\n%-34s\n%s" % (Name+' (spin):',FmtMat((a-b)/2,ind="    "))
            return r
         else:
            r = "%-34s\n%s" % ("%s (%s):" % (Name,Annotation),FmtMat(obj,ind="    "))
            return r
      else:
         r = "%-34s\n    %s" % (Name+':',obj)
   elif isinstance(obj, tuple):
      r = ("{:<34s}"+len(obj)*"{:14.8f}").format(Name,*obj)
   else:
      r = pResultFmtS % (Name,obj)
   if Annotation is not None:
      return r + "  %s" % Annotation
   else:
      return r
   pass

pInfoFmt = pResultFmt

def FmtHeading(s, isize):
   """print a subdivision marker between program subsections.
   isize == 0: very bold, isize == 1: very very bold."""
   N = 62
   if isize == 0:
      return "+ %s:" % s
   elif isize == 1:
      return (N/2) * "_ " + "\n\n  %s\n" % s.upper()
   elif isize > 1:
      return N * "_" + "\n\n  %s\n" % s + (N/2-1) * "_ " + "\n\n"

class FOutputLog(object):
   def __init__(self, OutputStream=stdout):
      self.out = OutputStream
      self.Level = 0
      self.Sections = []
      self.SectionLabels = [""]

      # can be used to enable/disable sections of output.
      self.PrintControl = {}

      self.SILENT = -1  # section is created, but not printed.
      self.ERROR = 0xff
      self.WARNING = 0x7f
      self.iInd = 0

   def __getitem__(self, p):
      # by default print everything. remember what was asked for.
      return self.PrintControl.setdefault(p, True)

   def _Emit(self, s):
      i = len(self.Sections)
      if i >= 10:
         i = 'a' + (i-10)
      Lead = "%1s%-20s %s"% (i,self.SectionLabels[-1],self.iInd*"  ")
      def EmitLine(s):
         print >> self.out, (Lead + s).rstrip()
      if "\n" not in s:
         EmitLine(s)
      else:
         Lines = s.split("\n")
         for Line in s.split("\n"):
            EmitLine(Line)
      pass

   def Ind(self):
      """increase indent-level (applied to all output)"""
      self.iInd += 1

   def UnInd(self):
      """decrease indent-level (applied to all output)"""
      self.iInd -= 1

   def Enter(self, Name, SectionDesc, iHeadingSize=0):
      self.Sections.append(Name)
      Separator = "%s/%s" if len(self.SectionLabels) != 1 else "%s:%s"
      Label = Separator % (self.SectionLabels[-1], Name)
      self.SectionLabels.append(Label)
      if iHeadingSize >= 0:
         self("{0}", FmtHeading(SectionDesc,iHeadingSize))
      self.Ind()

   def Leave(self, Name):
      self.NewLine()
      if self.Sections and self.Sections[-1] == Name:
         iDel = 1
      else:
         self(self.ERROR, "FOutputLog: leaving a non-top-level section. Current stack is:"
            " '{0}', trying to leave: '{1}'", "/".join(self.Sections), Name)
         # see if we have it somewhere further up the stack.
         S = self.Sections[::-1]
         if Name in S:
            iDel = S.index(Name) + 1
         else:
            # never entered -- nothing to do here.
            return
      self.Sections = self.Sections[:-iDel]
      self.SectionLabels = self.SectionLabels[:-iDel]
      self.UnInd()
      self.Flush()

   def Section(self, SectionName, SectionDesc, iHeadingSize=0):
      return FSectionLog(self, SectionName, SectionDesc, iHeadingSize)

   def NewLine(self):
      self("")

   def Flush(self):
      self.out.flush()

   def __call__(self, *p):
      """write output to the log, at the current logging level.
      Usage:
         L("wheee")
         L("{:10s} {:15.8f}", "something", 12.03)
      """
      if not p:
         return self.NewLine()
      assert(len(p) >= 1)
      Level = self.Level
      if isinstance(p[0], int):
         Level = p[0]
         p = p[1:]
      assert(len(p) >= 1)

      Fmt = p[0]

      if len(p) > 1:
         FmtArgs = p[1:]
         if hasattr(Fmt, '__call__'):
            # got passed a formatting function.
            Fmt = Fmt(*FmtArgs)
         else:
            # has to be a format string.
            assert(isinstance(Fmt,str))
            Fmt = Fmt.format(*FmtArgs)
      if Level == self.ERROR:
         Fmt = "\033[1;37;41mERROR:\033[0m " + Fmt
      if Level == self.WARNING:
         Fmt = "\033[1;37;43mWARNING:\033[0m " + Fmt
      self._Emit(Fmt)

class FSectionLog(object):
   def __init__(self, ParentLog, SectionName, SectionDesc, iHeadingSize):
      self.L = ParentLog
      self.Name = SectionName
      self.Desc = SectionDesc
      self.iHeadingSize = iHeadingSize
   def __enter__(self):
      self.L.Enter(self.Name, self.Desc, self.iHeadingSize)
   def __exit__(self, type, value, traceback):
      self.L.Leave(self.Name)


def FmtDateAndTime():
   import datetime
   import time
   #now = datetime.datetime.now()
   now = time.localtime()
   tzname = time.tzname[now.tm_isdst]
   # ^- now.tzname() doesn't work here.
   #    could also use time.daylight instead of tm_isdst.
   return "%04i-%02i-%02i %02i:%02i:%02i %s" % (tuple(now[:6]) + (tzname,))


def _TestOutputModule():
   import time
   import os
   L = FOutputLog()
   with L.Section("info", "Version & Runtime Environment", 1):
      L("{:30s} {}", "LatticeDmet", "v20121221")
      L(L.WARNING, "Sugar causes cavities. Take care.")
      L("{:30s} {}", "Current Time", FmtDateAndTime())
      (sysname, nodename, release, version, machine) = os.uname()
      L("{:30s} {} ({}, {} {})", "Running on", nodename, machine, sysname, release)

   L.Enter("lmf", "Lattice Mean-Field")
   L(pResultFmt, "Number of sites", 9000)
   L(pResultFmt, "Hartree-Fock energy", -230.73100019281212, "electronic")
   L(pResultFmt, "Hartree-Fock energy", -230.73100019281212, "total")
   L.Leave("lmf")

if __name__ == "__main__":
   _TestOutputModule()







