import numpy as np
import scipy.linalg as la
import settings as dmetSet
from os import remove, path

def Diag1eHamiltonian(Fock, nElec, ThrDeg, LastOrb = None, NumVirts = None, fOrbOcc = None):  
    nOrb = Fock.shape[0]
    assert(Fock.shape[0] == Fock.shape[1])
    assert(nElec <= nOrb)
    Ew, Orbs = la.eigh(Fock) 
    iEws = Ew.argsort()
    #print "Ew is:", Ew
    Ew = Ew[iEws]
    Orbs = Orbs[:,iEws]
    #print "difference is:", Fock-np.dot(Orbs*Ew,Orbs.T)
    #print "Ew is:", Ew
    if ( LastOrb is not None ):
        # rotate orbitals in degenerate subspace around the occupied edge
        # such that overlap with previous occupied space is maximized.
        iFirst = max(0, nOcc-40)
        iLast = min(N, nOcc+40)
        #ThrDeg = 1e-4
        while ( abs(ew[iLast-1] - ew[nOcc-1]) > ThrDeg ):
            iLast -= 1
        while ( abs(ew[iFirst] - ew[nOcc-1]) > ThrDeg ):
            iFirst += 1
        assert(iLast >= iFirst)
        if ( iLast - iFirst > 1 ):
            DegenSpace = orbs[:,iFirst:iLast]
            HalfProj = dot(transpose(LastOrb[:,:nOcc]),DegenSpace)
            OrderGen = mdot(transpose(HalfProj), HalfProj)
            ewg, evg = eigh(-OrderGen); ewg *= -1.
            if 0:
                print "!ew: %i--%i--%i    %s | %s | %s" % (iFirst,nOcc-1,iLast, ew[iFirst-3:iFirst], ew[iFirst:iLast], ew[iLast:iLast+3] )
                print "ewg:\n",ewg,"\nevg:\n",evg
             # rotate original orbitals into eigenvectors of order generator.
            Orbs[:,iFirst:iLast] = dot(DegenSpace, evg)
            ew[iFirst:iLast] = diag(mdot(evg.T, diag(ew[iFirst:iLast]), evg))
            # ^- should not do anything in principle. but the threshold above is rather large,
            #     and without it there may be severe convergence problems in some situations.
    Rdm = np.dot(Orbs[:,:nElec],Orbs[:,:nElec].T)
    Mu = (Ew[nElec-1] + Ew[nElec])/2.
    HlGap = Ew[nElec] - Ew[nElec-1]
    return Rdm, Mu, HlGap

def Diag1eGQNHamiltonian(Fock):  
    # gqn stands for good quantum number
    # Diagonalize 1e Hamiltonian which has good quantum number,
    # such as momentum k if theres' tranlational invarient
    Orbs = np.zeros_like(Fock)
    Ews = np.zeros((Fock.shape[0], Orbs.shape[2]))
    for gqn in range(Fock.shape[0]):
        ewi, evi = la.eigh(Fock[gqn,:,:])  
        Orbs[gqn,:,:] = evi
        Ews[gqn,:] = ewi
    return Ews,Orbs
    
def DealGQNOrbitals(Ews, Orbs, nElec, ThrDeg):  
    OccNumbers = np.zeros_like(Ews)
    Ews_ = Ews.flatten()
    iEws = Ews_.argsort()
    Lumo = Ews_[iEws[nElec]]
    Homo = Ews_[iEws[nElec-1]] if (nElec != 0) else Lumo
    Mu = .5 * (Homo + Lumo)
    HlGap = Lumo - Homo
    # first fill all orbitals below Mu-ThrDeg fully
    iOrbClo = Ews < Mu - ThrDeg
    iOrbOcc = Ews < Mu + ThrDeg
    nOrbClo = sum(iOrbClo.flatten())
    nOrbOcc = sum(iOrbOcc.flatten())
    nOrbAct = nOrbOcc - nOrbClo
    # now distribute the remaining electrons evenly amongst the degenerate
    # subspace
    nElecFull = nOrbClo
    nElecLeft = nElec - nElecFull
    assert(nElecLeft >= 0)
    if nOrbAct != 0:
        OccNumbers[iOrbOcc] = nElecLeft / nOrbAct
    OccNumbers[iOrbClo] = 1

    assert(np.all(OccNumbers >= 0.))
    # make the rdm: rdm[r,s] = \sum_i orb[r,i] occ-num[i] orb[s,i]
    #                                = \sum_i (orb[r,i] occ-num[i])**.5 * (orb[s,i] occ-num[i])**.5
    OrbsOcc = (Orbs.transpose((1,0,2)) * OccNumbers**.5).transpose((1,0,2))
    Rdm = np.zeros(Orbs.shape, OrbsOcc.dtype)
    for gqn in xrange(Rdm.shape[0]):
        C = OrbsOcc[gqn,:,:]    # <- nSites x nOrb
        Rdm[gqn,:,:] = np.dot(C, np.conj(C.T))
    return Rdm, Mu, HlGap


def MakeEmbeddingBasis(ImpSites, Rdm, ThrBathSvd):
    Rdm = Rdm.reshape((Rdm.shape[0]*Rdm.shape[1],) + Rdm.shape[2:]) 
    if (type(ImpSites) is not np.ndarray):
        ImpSites = np.array(ImpSites, int)
    nImp = len(ImpSites)
    EnvSites = [o for o in range(Rdm.shape[0]) if o not in ImpSites]
    EnvSites = np.array(EnvSites)
    EnvImpRdm = 1.*Rdm[EnvSites,:][:,ImpSites]
    S = np.dot(EnvImpRdm.T,EnvImpRdm)
    ew,ev = la.eigh(S)
    bathIndex = [list(ew).index(o) for o in ew if o > ThrBathSvd]
    nBath = len(bathIndex)
    bathBasis = np.dot(EnvImpRdm,ev[:,bathIndex]*ew[bathIndex]**(-0.5))
    print "the orthogonal of basis:\n", np.dot(bathBasis.T,bathBasis)
    embBasis = np.zeros((Rdm.shape[0],nImp+nBath))
    embBasis[ImpSites,:nImp] = np.eye(nImp)
    embBasis[EnvSites,nImp:] = bathBasis
    return embBasis


def MakeEmbeddingBasis1(ImpSites, Rdm, ThrBathSvd):
    if (type(ImpSites) is not np.ndarray):
        ImpSites = np.array(ImpSites, int)
    nImp = len(ImpSites)
    EnvSites = [o for o in range(Rdm.shape[0]) if o not in ImpSites]
    EnvSites = np.array(EnvSites)
    ImpRdm = 1.*Rdm[ImpSites,:][:,ImpSites]
    ew,ev = la.eigh(ImpRdm)
    bathIndex = [list(ew).index(o) for o in ew if o > ThrBathSvd]
    nBath = len(bathIndex)
    bathBasis = np.dot(Rdm[EnvSites,:][:,ImpSites],ev[:,bathIndex]*ew[bathIndex]**(-0.5))
    norm = np.zeros(nBath)
    for i in range(nBath):
        X = bathBasis[:,i]
        norm[i] = np.linalg.norm(X) 
    #print "the orthogonal of basis:\n", np.dot(bathBasis.T,bathBasis)
    bathBasis = bathBasis/norm#raise SystemExit
    print "the orthogonal of basis:\n", np.dot(bathBasis.T,bathBasis)
    embBasis = np.zeros((Rdm.shape[0],nImp+nBath))
    embBasis[ImpSites,:nImp] = np.eye(nImp)
    embBasis[EnvSites,nImp:] = bathBasis
    return embBasis

def ToEmb(M,EmbBasis):
    return np.dot(EmbBasis.T, np.dot(M,EmbBasis))

from textwrap import dedent

def WriteFile(FileName, Text):
    File = open(FileName, "w")
    File.write(Text)
    File.close()

def WriteDumpFile(FileName, Emb1e, n2eOrb, nElec, Ms2=0, U=None, Int2e=None, Uhf=False):
    if Uhf :
        nOrb = Emb1e.shape[0]/2
    else :
        nOrb = Emb1e.shape[0]
    Text = dedent("""
      &FCI NORB=%i,NELEC=%i,MS2=%i,
        ORBSYM=%s
        ISYM=1,%s
      &END
    """[1:-1])
    Text = Text % (nOrb, nElec, Ms2, nOrb*"1,",["","\n    IUHF=1,"][Uhf])
    Text = "\n".join(" " + o for o in Text.splitlines()) + "\n"
    IntLines = []
    assert(n2eOrb <= nOrb)
    LineFmt = "  %16.12f %i %i %i %i"
    def AppendSeparatorLine():
        # in UHF case, such lines separate between the different spin cases.
        IntLines.append(LineFmt % (0.0, 0, 0, 0, 0))
    def Append1eOp(Emb1e):
        for iOrb in range(nOrb):
            for jOrb in range(iOrb+1):
                tij = Emb1e[iOrb,jOrb]
                IntLines.append(LineFmt % (tij, iOrb+1,jOrb+1, 0, 0) )
    def Append2eOp(Int2e, U, iSymAB):
        #for iOrb in range(nImp):
        #    IntLines.append(LineFmt % (U, iOrb+1,iOrb+1,iOrb+1,iOrb+1) )
        if U is not None:
            for iOrb in range(n2eOrb):
                IntLines.append(LineFmt % (U, iOrb+1,iOrb+1,iOrb+1,iOrb+1) )
        elif Int2e is not None:
            assert(nOrb == Int2e.shape[0] == Int2e.shape[1] == Int2e.shape[2] == Int2e.shape[3])
            for i in range(nOrb):
                for j in range(nOrb):
                    for k in range(nOrb):
                        for l in range(nOrb):
                            f = Int2e[i,j,k,l]
                            i_ = max(i,j); j_ = min(i,j)
                            k_ = max(k,l); l_ = min(k,l)
                            ij = ((i_+1)*i_)/2+j_
                            kl = ((k_+1)*k_)/2+l_
                            if (i >= j and k >= l and (iSymAB == 0 or ij >= kl)):
                                IntLines.append(LineFmt % (f, i+1,j+1,k+1,l+1) )
    if ( not Uhf ):
        Append2eOp(Int2e,U,1)
        AppendSeparatorLine()
        Append1eOp(Emb1e)
    else:
        if ( Int2e is not None ):
            (Int2e_AA, Int2e_BB, Int2e_AB) = Int2e
        else:
            (Int2e_AA, Int2e_BB, Int2e_AB) = None, None, None
        Emb1eA = Emb1e[::2,::2]
        Emb1eB = Emb1e[1::2,1::2]
        Append2eOp(Int2e_AA,U,1)  # AA spin
        AppendSeparatorLine()
        Append2eOp(Int2e_BB,U,1)  # BB spin
        AppendSeparatorLine()
        Append2eOp(Int2e_AB,U,0)  # AB spin (does not have (ij|kl) = (kl|ij) symmetry)
        AppendSeparatorLine()
        Append1eOp(Emb1eA)
        AppendSeparatorLine()
        Append1eOp(Emb1eB)
        AppendSeparatorLine()
    IntLines.append(LineFmt % (0.0,0,0,0,0) ) # core energy.
    IntLines += [""] # last line: add a \n.
    Text = Text + "\n".join(IntLines)
    #print Text
    WriteFile(FileName, Text)

def OptionForHlMethod(HlMethod, FileNameRdm1,BasePath,IntType, nSys,  U=0,  nSysTrace=None,Int2e_=None,MakeRdm2=False,MakeRdm2S=False,MakeFock=False,MaxIt=2048,DiagBasis=None):
    if HlMethod == "Fci":     
        SubDim=16
        pSpaceDim = 200 # hm... doesn't play well toegether with spinproj.
        CiMethod = "FCI"
        if ( nSys < 8 ):
            pSpaceDim = 100
        if ( nSys >= 12 ):
            pSpaceDim = 1500
        ThrVar = 1e-12
        if DiagBasis is None:
            if ( U > 6. or (nSys > 4 and IntType != "RHF")):
                DiagBasis = "Input"
            else :
                DiagBasis = "CoreH"
                # ^- this can use the (barely working...) HubU optimization of fci.
        else :
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
        BaseCmd = " --subspace-dimension=%s --basis=%s --method='%s' --pspace=%i "\
            "--thr-var=%e --diis-block-size=%i --max-it=%i"\
            % (SubDim,DiagBasis,CiMethod,pSpaceDim,ThrVar,5000e3/2/SubDim,MaxIt)
        BaseCmd += " --save-rdm1='%s'" % FileNameRdm1
        if MakeFock:
            FileNameFock = path.join(BasePath, HlMethod+r"FOCK")
            BaseCmd += " --save-fock='%s'" % FileNameFock
        else:
            FileNameFock = None
        #BaseCmd += " --fci-vec='%s'" %FileNameFciVec
        if ( MakeRdm2 ):
            FileNameRdm2 = path.join(BasePath, HlMethod+r"2RDM")
            BaseCmd += " --save-rdm2='%s'" % FileNameRdm2
        if ( MakeRdm2S ):
            FileNameRdm2s = path.join(BasePath, HlMethod+r"2RDM-S")
            BaseCmd += " --save-rdm2s='%s'" % FileNameRdm2s
        if ( nSysTrace is not None ):
            BaseCmd = BaseCmd.replace("--thr-var", "--ptrace=%i --thr-var" % nSysTrace)
    else : 
        raise Exception("You have to implement the options for other HlMethods")
    return BaseCmd

def FindFciSolution(BasePath,nSys_,nSysTrace=None,CiMethod="FCI", IntType="RHF",
          U=0,Int2e_=None):
    FileNameInp = path.join(BasePath, r"FciINP")
    FileNameRdm = path.join(BasePath, r"FCIRdm")

    SubDim=16
    pSpaceDim = 200 # hm... doesn't play well toegether with spinproj.
    if ( U > 4. ):
        pSpaceDim = 400
    nSys = nSys_

    if ( nSys < 8 ):
        pSpaceDim = 100
    if ( nSys >= 12 ):
        pSpaceDim = 1500
    ThrVar = 1e-12
    #pSpaceDim = 0
    Executable = dmetSet.FciExecutable
    if dmetSet.param_fci["DiagBasis"] is None:
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
        else :
            WriteFile(FileNameBasis+"_A", "\n".join(MakeAsymOpLines(DiagBasis[ ::2, ::2])) + "\n")
            WriteFile(FileNameBasis+"_B", "\n".join(MakeAsymOpLines(DiagBasis[1::2,1::2])) + "\n")
            DiagBasis = "'!%s'" % FileNameBasis
    
    if dmetSet.param_fci["runOPENMP"] :
        BaseCmd = "OMP_NUM_THREADS=12 %d" % param_fci["nThreads"]
        BaseCmd += "%s --subspace-dimension=%s --basis=%s --method='%s' --pspace=%i "\
            "--thr-var=%e --diis-block-size=%i --max-it=%i"\
            % (Executable, SubDim,DiagBasis,CiMethod,pSpaceDim,ThrVar,5000e3/2/SubDim,dmetSet.param_fci["MaxIt"])
    else : 
        BaseCmd = "%s --subspace-dimension=%s --basis=%s --method='%s' --pspace=%i "\
            "--thr-var=%e --diis-block-size=%i --max-it=%i"\
            % (Executable, SubDim,DiagBasis,CiMethod,pSpaceDim,ThrVar,5000e3/2/SubDim,dmetSet.param_fci["MaxIt"])
    BaseCmd += " --save-rdm1='%s'" % FileNameRdm
    if dmetSet.param_fci["MakeFock"]:
        FileNameFock = path.join(BasePath, r"FOCK")
        BaseCmd += " --save-fock='%s'" % FileNameFock
    else:
        FileNameFock = None
    #BaseCmd += " --fci-vec='%s'" %FileNameFciVec
    if ( dmetSet.param_fci["MakeRdm2"] ):
        FileNameRdm2 = path.join(BasePath, r"2RDM")
        BaseCmd += " --save-rdm2='%s'" % FileNameRdm2
    if ( dmetSet.param_fci["MakeRdm2S"] ):
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
        Cmd = "%s '%s'" % (BaseCmd, FileNameInp)
        print "!%s" % Cmd
        Output = getoutput(Cmd)
        print Output
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
        print "fci Rdm1 is:\n", Rdm1
        if ( FileNameFock is not None ):
            FciFock = Read1Rdm(FileNameFock)
        #nOrb = h1.shape[0]
        if ( dmetSet.param_fci["MakeRdm2"] ):
            Rdm2 = Read1Rdm(FileNameRdm2).reshape((nOrb,nOrb,nOrb,nOrb))
        if ( dmetSet.param_fci["MakeRdm2S"] ):
            Rdm2S = Read1Rdm(FileNameRdm2s).reshape((nOrb,nOrb,nOrb,nOrb))
            #print Output
    except:
        Cleanup(BasePath)
        raise
    return  E, EpTrace, Rdm1, Cmd, Output
    #return FHlResults(Energy, EpTrace, nElec, Rdm1, Cmd, Output)

def readrdm(file, UHFB):
    with open(file, "r") as f:
        lines = f.readlines()
  
    nsites = int(lines[0])
    rdm = np.zeros((nsites, nsites))
  
    for line in lines[1:]:
        tokens = line.split(" ")
        rdm[int(tokens[0]), int(tokens[1])] = float(tokens[2])
  
    if UHFB:
        return [rdm[::2, ::2], rdm[1::2, 1::2]]
    else:
        return rdm

def FindDmrgSolution(BasePath, DmrgConfFile, IntType="RHF",p=True):
    #cmd = ["mpiexec -np %d" % nproc]
    #if param_block["bind"]:
    #  cmd.append("--bind-to-socket")
    #cmd.append(param_block["exec"])
    #cmd.append("dmrg.conf")
    #cmd = " ".join(cmd)
    # call block
    FileNameInp = path.join(BasePath, r"DmrgINP")
    DmrgOutFile = path.join(BasePath, r"DmrgOUT")
    if dmetSet.param_block["runMPI"] :
        BaseCmd = "mpiexec -np %d" %dmetSet.param_block["nproc"]
        if dmetSet.param_block["bind"]:
            BaseCmd += " --bind-to-socket" 
        BaseCmd += dmetSet.DmrgExecutable
    else : 
        BaseCmd = dmetSet.DmrgExecutable
    Cmd = "%s %s > %s " %(BaseCmd, DmrgConfFile, DmrgOutFile)
    print Cmd
    Output = getoutput(Cmd)
    Output = getoutput("grep \"Sweep Energy\" %s"%(DmrgOutFile))
    if p:
        print Output
        print
    # read energy
    #file_e = open("dmrg.e", "rb")
    file_e = open(path.join(BasePath,r"dmrg.e"),"rb")
    energy = struct.unpack('d', file_e.read(8))[0]
    EpTrace = 0.0
    # read rdm and kappa
    #if IntType == "RHF":
    Rdm1 = readrdm("spatial_onepdm.0.0", False)
    #if IntType == "UHF":
    #    rdm = readrdm("onepdm.0.0.txt", True)
    #else:
    #    rdm = readrdm("spatial_onepdm.0.0.txt", False) / 2
    #    kappa = readrdm("spatial_pairmat.0.0.txt", False)
    #if not UHFB:
    #    kappa = (kappa + kappa.T) / 2
    #os.chdir(cwd)
    #cwd = os.getcwd()
    #os.chdir(path)
    Cleanup(BasePath)
    return  E, EpTrace, Rdm1, Cmd, Output

def WriteDmrgConfigFile(BasePath, nsites,nelec, Ms2=0):
    configfile = []
    blockM = dmetSet.param_block["M"]
    if nelec == -1:
        configfile.append("nelec %d\nspin 0\nhf_occ integral" % (nsites*2))
        configfile.append("bogoliubov" % (nsites*2))
    else:
        #configfile.append("nelec %d\nspin %d\n hf_occ integral" % (nelec, Ms2))
        configfile.append("nelec %d\nspin %d" %(nelec, Ms2))

    configfile.append(block_schedule(blockM))

    configfile.append("orbitals "+BasePath+"/DmrgINP")
    configfile.append("prefix  " +BasePath)
    #configfile.append("nonspinadapted\nonepdm")
    configfile.append("onepdm")
    #configfile.append("noreorder")
    configfile.append("outputlevel -1")
    configtxt = "\n".join(configfile) + "\n"
    configFilename = BasePath+"/dmrg.conf" 
    WriteFile(configFilename, configtxt)
    return configFilename
  

def block_schedule(M):
    lines = ["schedule"]
    if M == 200:
        lines.append("0 50  1e-5 1e-4")
        lines.append("4 100 1e-5 1e-4")
        lines.append("8 200 1e-6 1e-5")
        lines.append("10 200 1e-7 0")
        lines.append("end")
        lines.append("")
        lines.append("twodot_to_onedot 12")
        lines.append("sweep_tol 1e-6")
    elif M == 300:
        lines.append("0 50  1e-5 1e-4")
        lines.append("4 100 1e-5 1e-4")
        lines.append("8 250 2e-6 1e-5")
        lines.append("12 300 1e-6 1e-5")
        lines.append("14 300 1e-7 0")
        lines.append("end")
        lines.append("")
        lines.append("twodot_to_onedot 16")
        lines.append("sweep_tol 1e-6")
    elif M == 400:
        lines.append("0 50  1e-5 1e-4")
        lines.append("4 100 1e-5 1e-4")
        lines.append("8 250 2e-6 1e-5")
        lines.append("12 400 1e-6 1e-6")
        lines.append("14 400 1e-7 0")
        lines.append("end")
        lines.append("")
        lines.append("twodot_to_onedot 16")
        lines.append("sweep_tol 1e-6")
    else:
        print "Block schedule not specified!"
        assert()

    lines.append("maxiter 30")
    return "\n".join(lines)

def ReadFile(FileName):
    File = open(FileName, "r")
    Text = File.read()
    File.close()
    return Text

def Read1Rdm(FileName):
    Text = ReadFile(FileName)
    Lines = Text.splitlines()
    nRows,x,nCols = Lines[0].split()[-3:]
    nRows = int(nRows)
    nCols = int(nCols)
    Numbers = map(float,(" ".join(Lines[1:])).split())
    return np.array(Numbers).reshape(nRows,nCols)

def FindOutput(Output,s):
    Lines= Output.splitlines()
    i = len(Lines) - 1
    while ("!%s STATE 1 %s" % ("FCI",s)) not in Lines[i]:
        i -= 1
    E = float(Lines[i].split()[-1])
    return E

from shutil import rmtree
def Cleanup(BasePath):
    rmtree(BasePath)
    pass

import os
from commands import getoutput
if __name__ == '__main__':
    nOcc = 3
    nOrb = 8
    nk = 5
    #Fock = np.zeros((nk,nOrb,nOrb))
    #for i in range(nk):
    #    Fock[i,:,:] = np.random.random((nOrb,nOrb))
    #    Fock[i,:,:] = Fock[i,:,:] + Fock[i,:,:].T
    #Ews, Orbs = Diag1eGQNHamiltonian(Fock)  
    #Rdm, Mu, Gap =  DealGQNOrbitals(Ews, Orbs, nOcc, 1.0e-3)  
    #print Rdm.shape

    Fock = np.random.random((nOrb,nOrb))
    Fock = Fock + Fock.T
    Rdm, Mu, Gap = Diag1eHamiltonian(Fock, nOcc, 1.0e-6) 
    ew,ev = la.eigh(Fock)
    Rdm = np.dot(ev[:,:nOcc],ev[:,:nOcc].T)
    ImpSites = [0,3]
    ThrBathSvd = 1.0e-8
    EmbBasis = MakeEmbeddingBasis(ImpSites, Rdm, ThrBathSvd)
    EmbBasis = MakeEmbeddingBasis1(ImpSites, Rdm, ThrBathSvd)
    EmbFock = ToEmb(Fock,EmbBasis)
    EmbRdm = ToEmb(Rdm,EmbBasis)
    nElec = int(np.trace(EmbRdm))
    print "nElec is", nElec
    print "EmbFock is:\n", EmbFock
    EmbRdm = Diag1eHamiltonian(EmbFock, nElec, 1.0e-6) 
    print "EmbRdm is:\n", 2*EmbRdm
    #WriteFciInput("./FciInput", EmbFock, 2*nElec, U=0.0 )
    #Output = getoutput("../fci.20121221/fci FciInput --save-rdm1 rdm1 FciInput")
    #print Output
    #Rdm1 = Read1Rdm("rdm1")
    #print "Rdm1 is:\n", Rdm1

