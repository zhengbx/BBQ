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

from numpy import *
import numpy as np
from numpy.linalg import *
import scipy
import scipy.optimize


def PrintHeading(s, isize):
   """print a subdivision marker between program subsections.
   isize == 0: very bold, isize == 1: very very bold."""
   N = 62
   if ( isize == 0 ):
      print (N/2) * "_ " + "\n\n %s\n" % s
   else:
      print N * "_" + "\n\n %s\n" % s + (N/2) * "_ " + "\n\n"

#def PrintMatrix(Text, M):
   #import settings as g
   #ESC = ExtractSpinComp
   #if ( len(M.shape) == 1 ):
      #if ( g.WF_TYPE == "RHF" ):
         #print " %s:\n %s" % (Text,M)
      #else:
         ##print "%s:\nA%s\nB%s\n" % (Text, ESC(M,0), ESC(M,1))
         #print " %s:\n C%s\n S%s\n" % (Text, ESC(M,0)+ESC(M,1), ESC(M,0)-ESC(M,1))
   #elif ( len(M.shape) == 2 ):
      #if ( g.WF_TYPE == "RHF" ):
         #print "%s:\n%s" % (Text,M)
      #else:
         #print "%s [charge]:\n%s\n%s [spin]:\n%s\n"\
            #% (Text, ESC(M,0)+ESC(M,1), Text, ESC(M,0)-ESC(M,1))
   #else:
      #assert(0)

def MakeSmh(S):
   """calculate S^{-1/2}."""
   ew,ev = eigh(-S); ew *= -1
   if ( ew[-1] < 1e-8 ):
      print ew
      raise Exception("S^{-1/2} ill-conditioned. Smallest eigenvalue: %8.2e" % ew[-1])
   evx = ev * (ew**-.25)
   return dot(evx, evx.T)

def mdot(*args):
   """chained matrix product."""
   r = args[0]
   for a in args[1:]:
      r = dot(r,a)
   return r

def dot2(X,Y):
   assert(X.shape == Y.shape)
   return dot(X.flatten(), Y.flatten())

def ReadFile(FileName):
   File = open(FileName, "r")
   Text = File.read()
   File.close()
   return Text

def WriteFile(FileName, Text):
   File = open(FileName, "w")
   File.write(Text)
   File.close()

def Read1Rdm(FileName):
   Text = ReadFile(FileName)
   Lines = Text.splitlines()
   nRows,x,nCols = Lines[0].split()[-3:]
   nRows = int(nRows)
   nCols = int(nCols)
   Numbers = map(float,(" ".join(Lines[1:nRows+1])).split())
   return array(Numbers).reshape(nRows,nCols)

def is_square(h):
    return len(h[0,:]) == len(h[:,0])

def is_hermitian(X):
   return np.allclose(X - conj(X.T), 0.0)

def to_square(mat, n = None):
    if ( n is None ):
        n = int(round(((1+8*len(mat))**.5 - 1)/2))
    square_mat=zeros([n,n],mat.dtype)
    ptr=0
    for i in xrange(n):
        for j in xrange(i+1):
            square_mat[i,j]=mat[ptr]
            square_mat[j,i]=square_mat[i,j]
            ptr+=1
    return square_mat

def to_triangle(mat):
    assert(mat.shape[0] == mat.shape[1])
    n = mat.shape[0]
    r = zeros(((n+1)*n)/2,mat.dtype)
    k = 0
    for i in xrange(n):
        for j in xrange(i+1):
            r[k] = mat[i,j]
            k += 1
    assert(sum((to_square(r,n) - mat)**2) < 1e-10)
    return r

def InvertPermutation(P):
   IP = len(P) * [0]
   for (i,p) in enumerate(P):
      IP[p] = i
   if type(P) is np.ndarray:
      return array(IP)
   else:
      return IP

IsSquare = is_square
ToSquare = to_square
ToTriangle = to_triangle

def ElementNameDummy():
   ElementNames = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn".split()
   ElementNumbers = dict([(o,i+1) for (i,o) in enumerate(ElementNames)])
   return ElementNames, ElementNumbers
ElementNames, ElementNumbers = ElementNameDummy()

def ExtractSpinComp(h,iComp, nSpinComp=2):
    # return alpha or beta components
    if ( h is None ):
        return None
    if ( nSpinComp == 1 ):
        return h
    if (len(h.shape) == 1):
        assert(h.shape[0] % 2 == 0)
        return h[iComp::nSpinComp]
    assert(len(h.shape) == 2)
    assert(h.shape[0] % 2 == 0 and h.shape[1] % 2 == 0)
    assert(iComp < nSpinComp)
    return h[iComp::2,iComp::2]

def CombineSpinComps(hAlpha, hBeta, nSpinComp=2):
    if ( hAlpha is None ):
        return None
    if ( nSpinComp == 1 ):
        return hAlpha
    assert(hAlpha.shape == hBeta.shape)
    if ( len(hAlpha.shape) == 2 ):
        N,M = hAlpha.shape
        h = zeros((2*N,2*M),hAlpha.dtype)
        h[::2,::2] = hAlpha
        h[1::2,1::2] = hBeta
        return h
    elif ( len(hAlpha.shape) == 1 ):
        N = hAlpha.shape[0]
        h = zeros(2*N,hAlpha.dtype)
        h[::2] = hAlpha
        h[1::2] = hBeta
        return h
    else:
        assert(0)


def resmin(fn, x0, Algo="Hard"):
    from numpy.linalg import solve, lstsq
    f0 = fn(x0)
    assert(len(f0.shape) == 1)
    nr = f0.shape[0]    # num residuals
    nx = len(x0)        # num variables
    def MakeGradMatrix(x):
        g = zeros((nr,nx))
        step = 1e-7
        for ix in range(nx):
            dx = zeros((nx))
            dx[ix] = step
            g[:,ix] = (.5/step)*(fn(x0 + dx) - fn(x0 - dx))
        return g
    x = x0

    p = False
    if 1:
        #Diis = FDiisContext(5)
        xval = dot(x,x)**.5
        step = 3e-3
        for it in range(200):
        #for it in range(20):
            r = fn(x)
            Err = dot(r,r)
            #if ( it >= 1 and False ):
                #x,r,diis_c0 = Diis.Apply(x,r)
                #r = fn(x)
            if ( Err < 1e-15 and it != 0 ):
                if p: print "       |%3i  %.3e" % (it,Err)
                break
            if Algo == "Hard" or it == 0:
               g = MakeGradMatrix(x)
            else:
               # bfgs update
               y = r - LastR
               s = dx
               gs = dot(g,s)
               g += outer(y,y)/dot(y,s) - outer(gs,gs)/dot(gs,s)

            if Algo == "Perturbative":
               (dx,fitresid,rank,sigma) = lstsq(g, r)
               return -dx

            def GetDir():
                #gsq = dot(transpose(g),g)
                #ew,ev = eigh(gsq)
                #print "gsq/ew:", ew
                #return ev[:,-1]

                #if nr != nx:
                if 0:
                    # get an update for x.
                    # we solve r == g * dx for dx.
                    (dx,fitresid,rank,sigma) = lstsq(g, r)
                else:
                    # get direction -- with error-dependent damping for
                    # step length restriction
                    N = len(r)
                    g2 = zeros((nr+nx,nx))
                    g2[:nr,:nx] = g
                    g2[nr:,:nx] = .1*Err**.5 * eye(nx)
                    #g2[nr:,:nx] = 5.*eye(N)/len(x)**2
                    r2 = zeros((nr+nx))
                    r2[:nr] = r
                    (dx2,fitresid,rank,sigma) = lstsq(g2, r2)
                    dx = dx2[:nx]
                return dx
            dx = GetDir()

            #print "g:\n",g
            #print "r:\n",r
            #print "dx:\n",dx

            def LineSearchFn(step):
                err1 = fn(x - step * dx)
                return dot(err1,err1)

            def FindStep():
                if Algo == "Simple":
                    return 1.
                elif True:
                    grid = list(arange(-0.0,2.00001,0.2))
                    #grid = list(arange(-0.0,2.00001,0.02))
                    #grid += list(arange(-0.0,0.02,0.0002))
                    val = [LineSearchFn(step) for step in grid]
                    s = grid[argmin(val)]
                    #return s
                    if ( abs(s) > 1e-4 ):
                        return s
                    else:
                        return scipy.optimize.fmin(LineSearchFn, array([0.001]),disp=0,xtol=1e-10)
                else:
                    step1 = scipy.optimize.fmin(LineSearchFn, array([1.000]),disp=0,xtol=1e-10)
                    step2 = scipy.optimize.fmin(LineSearchFn, array([0.001]),disp=0,xtol=1e-10)
                    if ( LineSearchFn(step1) <= LineSearchFn(step2) ):
                        return step1
                    else:
                        return step2
            step = FindStep()
            #step = 1.0

            #while ( abs(step) < 1e-4 ):
            for i in range(0):
                x -= step * dx
                ndx = dx/norm(dx)
                g = g - outer(ndx,dot(ndx,g))
                dx = GetDir()
                step = FindStep()
                print "           new direction:  step = %.3e" % step

            if p: print "       |%3i  %.3e   %.3e  %+.2e  %s.." % (it,dot(r,r),norm(dx),step, x[:5])

            if ( abs(step)*norm(dx) < 1e-10 ):
                # take the first local minimum
                break

            dx *= step
            x -= dx
            LastR = r
            #x,r,diis_c0 = Diis.Apply(x,fn(x))
        #if ( norm(r) > 1e-10 ):
           #print "WARNING: resmin failed to converge. Final gradient: %8.2e" % norm(r)
        return x
