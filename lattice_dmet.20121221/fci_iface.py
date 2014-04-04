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

# -*- coding: utf-8 -*-
from scipy.linalg import *
from scipy import *

from textwrap import dedent
from commands import getoutput
from os import remove, path
from tempfile import mkdtemp
from shutil import rmtree

import settings as g
from helpers import *


def WriteFciInput(FileName,h,U,nSys,nElec,Ms2=0,h1b=None,Int2e=None,nOrb2e=None):
    nOrb = len(h[:,0])

    # h1b, if provided, is a separate operator for beta spin one-electron
    # integrals. In that case we carry out the FCI in spin-orbital basis
    # (instead of spatial orbital basis). File format looks a bit different
    # then.
    Uhf = h1b is not None
    Text = dedent("""
      &FCI NORB=%i,NELEC=%i,MS2=%i,
       ORBSYM=%s
       ISYM=1,%s
      &END
    """[1:-1])
    Text = Text % (nOrb, nElec, Ms2, nOrb*"1,",["","\n   IUHF=1,"][Uhf])
    Text = "\n".join(" " + o for o in Text.splitlines()) + "\n"
    #WriteFile(FileName,Text)
    IntLines = []
    assert(nSys <= nOrb)
    #for iOrb in range(nSys):
        ## that's [1/2] e^ii_ii = nalpha*nbeta, if I'm not mistaken.
        #IntLines.append("  %16.12f %i %i %i %i" % (U, iOrb+1,iOrb+1,iOrb+1,iOrb+1) )

    LineFmt = "  %16.12f %i %i %i %i"
    def AppendSeparatorLine():
        # in UHF case, such lines separate between the different spin cases.
        IntLines.append(LineFmt % (0.0, 0, 0, 0, 0))
    def Append2eOp(Int2e,U, iSymAB):
        # iSymAB: if set, only unique integrals under the condition
        # (ij|kl) = (kl|ij) are emitted.
        if ( Int2e is None ):
            # write standard hubbard integrals.
            nOrb2e_ = nOrb2e
            if ( nOrb2e_ is None ):
                nOrb2e_ = nSys
            for iOrb in range(nOrb2e_):
                # that's [1/2] e^ii_ii = nalpha*nbeta, if I'm not mistaken.
                IntLines.append(LineFmt % (U, iOrb+1,iOrb+1,iOrb+1,iOrb+1) )
        else:
            # write input integrals.
            assert(nOrb == Int2e.shape[0] == Int2e.shape[1] == Int2e.shape[2] == Int2e.shape[3])
            for i in range(nOrb):
                for j in range(nOrb):
                    for k in range(nOrb):
                        for l in range(nOrb):
                            f = Int2e[i,j,k,l]
                            #if ( abs(f) > 1e-12 ):
                            # ^- FIXME: put this in and everything blows up and
                            # gives random results. why?
                            i_ = max(i,j); j_ = min(i,j)
                            k_ = max(k,l); l_ = min(k,l)
                            ij = ((i_+1)*i_)/2+j_
                            kl = ((k_+1)*k_)/2+l_
                            if (i >= j and k >= l and (iSymAB == 0 or ij >= kl)):
                                IntLines.append(LineFmt % (f, i+1,j+1,k+1,l+1) )
                            #if ( i >= j and k >= l and (iSymAB == 0 or ((i+1)*i)/2+j >= ((k+1)*k)/2+l)):
                                #IntLines.append(LineFmt % (f, i+1,j+1,k+1,l+1) )

    def Append1eOp(h):
        for iOrb in range(nOrb):
            for jOrb in range(iOrb+1):
                tij = h[iOrb,jOrb]
                #if ( abs(tij) > 1e-14 ):
                if True:
                    IntLines.append(LineFmt % (tij, iOrb+1,jOrb+1, 0, 0) )
    if ( not Uhf ):
        Append2eOp(Int2e,U,1)
        Append1eOp(h)
    else:
        if ( Int2e is not None ):
            (Int2e_AA, Int2e_BB, Int2e_AB) = Int2e
        else:
            (Int2e_AA, Int2e_BB, Int2e_AB) = None, None, None
        Append2eOp(Int2e_AA,U,1)  # AA spin
        AppendSeparatorLine()
        Append2eOp(Int2e_BB,U,1)  # BB spin
        AppendSeparatorLine()
        Append2eOp(Int2e_AB,U,0)  # AB spin (does not have (ij|kl) = (kl|ij) symmetry)
        AppendSeparatorLine()
        Append1eOp(h)
        AppendSeparatorLine()
        Append1eOp(h1b)
        AppendSeparatorLine()


    IntLines.append(LineFmt % (0.0,0,0,0,0) ) # core energy.
    IntLines += [""] # last line: add a \n.
    Text = Text + "\n".join(IntLines)
    #print Text
    WriteFile(FileName, Text)

def Read1Rdm(FileName):
    Text = ReadFile(FileName)
    Lines = Text.splitlines()
    nRows,x,nCols = Lines[0].split()[-3:]
    nRows = int(nRows)
    nCols = int(nCols)
    Numbers = map(float,(" ".join(Lines[1:])).split())
    return array(Numbers).reshape(nRows,nCols)


class DivergenceError(Exception):
    pass

#FileNameFciVec = path.join(g.TmpDir, r"FCIVEC")
#FileNameInp1 = None

def FindFciSolution(h1_,U,nSys_,nElec,CiMethod,IntType="RHF",Ms2=0,
        nSysTrace=None,Int2e_=None,MakeRdm2=False,MakeRdm2S=False,MakeFock=False,MaxIt=2048,
        DiagBasis=None):
    #for x in "h1_ DiagBasis Int2e_".split():
        #print "!xa2 %s.shape = %s" % (x, locals()[x].shape)
    #raise SystemExit
    BasePath = mkdtemp(prefix="fci", dir=g.TmpDir)
    def Cleanup():
       rmtree(BasePath)
       pass
    FileNameInp = path.join(BasePath, r"FCIINP")
    FileNameRdm = path.join(BasePath, r"1RDM")
    global FileNameInp1
    FileNameInp1 = FileNameInp
    #print "!setting FileNameInp1 = %s" % FileNameInp1

    if ( IntType == "RHF" ):
        h1 = h1_
        h1b = None
        Int2e = Int2e_
        nSys = nSys_
    else:
        assert(IntType == "UHF")
        # 2-RDMs in UHF mode currently not supported by fci program
        #   (but can be coded if required). Would also need some adjustements
        #   below for re-combining the spatial spin-parts into the spin-orb
        #   basis.
        assert((MakeRdm2 is None or MakeRdm2 == False))
        assert((MakeRdm2S is None or MakeRdm2S == False))
        # nSys_ is given in terms of spin-orbitals. Convert to spatial.
        assert(nSys_ % 2 == 0)
        nSys = nSys_/2

        # extract alpha and beta components of 1e Hamiltonian.
        h1  = h1_[0::2,0::2]
        h1b = h1_[1::2,1::2]
        # cross components connecting alpha- and beta orbitals are
        # not supported in this FCI program.
        assert( norm(h1_[0::2,1::2]) < 1e-10 )
        if Int2e_ is None:
            Int2e = None
        else:
            # get tuple for AA, BB, and AB components of 2e interaction.
            Int2e = (Int2e_[0::2,0::2,0::2,0::2],
                     Int2e_[1::2,1::2,1::2,1::2],
                     Int2e_[0::2,0::2,1::2,1::2])

    SubDim=16
    pSpaceDim = 200 # hm... doesn't play well toegether with spinproj.
    if ( U > 4. ):
        pSpaceDim = 400

    if ( h1.shape[0] < 8 ):
        pSpaceDim = 100
    if ( h1.shape[0] >= 12 ):
        pSpaceDim = 1500
    ThrVar = 1e-12
    #pSpaceDim = 0
    Executable = g.FciExecutable
    if ( CiMethod == "DMRG" ):
        if ( h1.shape[0] <= 2 or nElec <= 2 ):
            # DMRG is buggy for less than two orbitals or two sites.
            CiMethod = "FCI"
        else:
            Executable = DmrgExecutable
            #ThrVar *= 1e2 # dmrg has serious troubles with hard thresholds.
            #ThrVar *= 1e-2
    if DiagBasis is None:
        if ( U > 6. or (nSys > 4 and IntType != "RHF")):
            DiagBasis = "Input"
        else:
            DiagBasis = "CoreH"
            # ^- this can use the (barely working...) HubU optimization of fci.
    else:
        if type(DiagBasis) is not str:
            # it's an array. write it to disk
            # (used for ability to perform limited CI in canonical HF basis)
            BasisLines = []
            #assert(allclose(dot(DiagBasis.T,DiagBasis) - eye(DiagBasis.shape[0]),0.))
            def MakeAsymOpLines(Op):
                assert(allclose(dot(Op.T,Op) - eye(Op.shape[0]),0.))
                Lines = []
                for j in range(Op.shape[1]):
                    for i in range(Op.shape[0]):
                        Lines.append("%6i %6i %24.15e" % (1+i,1+j,Op[i,j]))
                return Lines
            FileNameBasis = path.join(BasePath, r"1DIAG_BASIS")
            if IntType == "RHF":
                WriteFile(FileNameBasis, "\n".join(MakeAsymOpLines(DiagBasis)) + "\n")
            else:
                WriteFile(FileNameBasis+"_A", "\n".join(MakeAsymOpLines(DiagBasis[ ::2, ::2])) + "\n")
                WriteFile(FileNameBasis+"_B", "\n".join(MakeAsymOpLines(DiagBasis[1::2,1::2])) + "\n")
            DiagBasis = "'!%s'" % FileNameBasis

    BaseCmd = "%s --subspace-dimension=%s --basis=%s --method='%s' --pspace=%i "\
        "--thr-var=%e --diis-block-size=%i --max-it=%i"\
        % (Executable, SubDim,DiagBasis,CiMethod,pSpaceDim,ThrVar,5000e3/2/SubDim,MaxIt)
    BaseCmd += " --save-rdm1='%s'" % FileNameRdm
    if MakeFock:
        FileNameFock = path.join(BasePath, r"FOCK")
        BaseCmd += " --save-fock='%s'" % FileNameFock
    else:
        FileNameFock = None
    #BaseCmd += " --fci-vec='%s'" %FileNameFciVec
    if ( MakeRdm2 ):
        FileNameRdm2 = path.join(BasePath, r"2RDM")
        BaseCmd += " --save-rdm2='%s'" % FileNameRdm2
    if ( MakeRdm2S ):
        FileNameRdm2s = path.join(BasePath, r"2RDM-S")
        BaseCmd += " --save-rdm2s='%s'" % FileNameRdm2s

    #BaseCmd += " --thr-print-c=0."

    if ( nSysTrace is not None ):
        BaseCmd = BaseCmd.replace("--thr-var", "--ptrace=%i --thr-var" % nSysTrace)

    def FindOutput(Output,s):
        Lines= Output.splitlines()
        i = len(Lines) - 1
        while ("!%s STATE 1 %s" % (CiMethod,s)) not in Lines[i]:
            i -= 1
        E = float(Lines[i].split()[-1])
        return E
    try:
        WriteFciInput(FileNameInp,h1,U,nSys,nElec,h1b=h1b,Ms2=Ms2,Int2e=Int2e)
        Cmd = "%s '%s'" % (BaseCmd, FileNameInp)
        #print "!%s" % Cmd
        Output = getoutput(Cmd)
        try:
            E = FindOutput(Output,"ENERGY")
            EpTrace = None
            if nSysTrace is not None:
                EpTrace = FindOutput(Output,"pTraceSys")
        except IndexError,e:
            print "Calculation failed: Output = ", Output
            print "Command was:"
            print "!%s" % Cmd

        # read 1-RDM.
        if ( IntType == "RHF" ):
            Rdm1 = Read1Rdm(FileNameRdm)
        else:
            Rdm1a = Read1Rdm(FileNameRdm + ".A")
            Rdm1b = Read1Rdm(FileNameRdm + ".B")
            Rdm1 = CombineSpinComps(Rdm1a, Rdm1b)
        HL = {\
            "Energy": E,
            "nElec": nElec,
            "Rdm1": Rdm1,
            "EpTrace": EpTrace,
            "Output": Output,
            "Cmd": Cmd,
        }
        if ( FileNameFock is not None ):
            HL["FciFock"] = Read1Rdm(FileNameFock)
        nOrb = h1.shape[0]
        if ( MakeRdm2 ):
            HL["Rdm2"] = Read1Rdm(FileNameRdm2).reshape((nOrb,nOrb,nOrb,nOrb))
        if ( MakeRdm2S ):
            HL["Rdm2S"] = Read1Rdm(FileNameRdm2s).reshape((nOrb,nOrb,nOrb,nOrb))
            #print Output
    except:
        Cleanup()
        raise

    return HL


# kate: indent-mode normal; indent-width 4; tab-width 4;
