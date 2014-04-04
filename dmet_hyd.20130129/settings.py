# This program is free software. It comes without any warranty, to the extent
# permitted by applicable law. You may use it, redistribute it and/or modify
# it, in whole or in part, provided that you do so at your own risk and do not
# hold the developers or copyright holders liable for any claim, damages, or
# other liabilities arising in connection with the software.
# 
# Developed by Gerald Knizia and Garnet K.-L. Chan, 2012;
# (c) Princeton University, 2012

from os import path
MakeIntegralsExecutable = "/home/cgk/dev/wmme/wmme"
FciExecutable = "/home/cgk/dev/fci/fci"

TmpDir = None
#TmpDir = "/dev/shm/"
# ^- if not None, use this as prefix for temporary files. Otherwise
#    use whatever mkdtmp() generates (likely related to TMPDIR environment
#    variable)
#BasisLibDir = "/home/cgk/dev/ct8k/"
BasisLibDir = path.join(path.dirname(MakeIntegralsExecutable),"bases")

from numpy import set_printoptions, nan
set_printoptions(precision=8,linewidth=10060,suppress=True,threshold=nan)


pResultFmt = " %-32s%18.12f"
pResultFmtAnnoted = " %-32s%18.12f  (%s)"
pResultFmtS = " %-32s%18s"
pResultFmtI = " %-32s%5i"
pResultFmtIAnnoted = " %-32s%5i  (%s)"


THRORB = 1e-10      # threshold for orbital gradient in RHF
THRDEN = 1e-8       # threshold for change in energy in RHF
MAXIT_HF = 10002     # maximum num. iterations for RHF

THRDVC = 1e-6       # threshold for change in full system correlation potential
MAXIT_VC = 128      # maximum num. iterations for correlation potential

IPRINT = 0          # 3: print full sys & sub.sys calculations
                    # 1: print report on sub.sys for each vcorr iteration
                    # 0: print vcorr summary and sub.sys info only after final iteration.

VLOC_FIT_TYPE = ["ImpRdm", "FullRdm", "EnvRdm", "EnvAndCoupling", "Coupling", "ImpRdmAndEnvRdm"][0]

                    # controls whether breaking of spin-symmetry is allowed.
WF_TYPE = ["RHF","UHF"][0]

VLOC_TYPE = ["Local","ImpAndEnv","Diagonal"][0]

USE_INT2E = True    # decides between interacting and non-interacting embedding.

ABSORB_HF_CORE = False # if true, start with a RHF calculation on the full system and
                       # remove core orbitals and electrons from all further considerations.
                       # Does not work with UHF.


ToAng2006 = 0.529177249   # 2006 value.
ToAng =     0.5291772108  # molpro default.
ToAngCp2k = 0.529177211   # CP2K default

# occupation number of an occupied orbital
if WF_TYPE == "RHF":
   ORB_OCC = 2.
else:
   assert(WF_TYPE == "UHF")
   ORB_OCC = 1.


