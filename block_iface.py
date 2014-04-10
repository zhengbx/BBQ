import numpy as np
import itertools as it
from shutil import rmtree
from commands import getoutput
from tempfile import mkdtemp
import os
import struct

from helpers import WriteFile, ReadFile
from settings import TmpDir, param_block

def empty_line():
  return "%20.16f%4d%4d%4d%4d" % (0., 0, 0, 0, 0)

def insert_ccdd(lines, v, nsites, seperation = True):
  # (ij||kl) convention indices are i k j l
  for idx1, (i, k) in enumerate(it.product(range(nsites), repeat = 2)):
    for idx2, (j, l) in enumerate(it.product(range(nsites), repeat = 2)):
      if v.has_key((i,j,l,k)):
        lines.append("%20.16f%4d%4d%4d%4d" % (v[(i,j,l,k)], i+1,k+1,j+1,l+1))
        del v[(i,j,l,k)]
  assert(len(v) == 0)
  if seperation:
    lines.append(empty_line())

def insert_cccd(lines, v, nsites):
  # plain index i j k l
  for i in range(nsites):
    for j in range(i):
      for k, l in it.product(range(nsites), repeat = 2):
        if v.has_key((i,j,k,l)):
          lines.append("%20.16f%4d%4d%4d%4d" % (v[(i,j,k,l)], i+1,j+1,k+1,l+1))
          del v[(i,j,k,l)]
  assert(len(v) == 0)
  lines.append(empty_line())

def insert_cccc(lines, v, nsites, UHFB):
  # index order i j l k, where i>j, l>k
  for i in range(nsites):
    for j in range(i):
      if UHFB:
        for l in range(nsites):
          for k in range(l):
            if v.has_key((i,j,k,l)):
              lines.append("%20.16f%4d%4d%4d%4d" % (v[(i,j,k,l)], i+1,j+1,l+1,k+1))
              del v[(i,j,k,l)]
      else:
        for l in range(i+1):
          for k in range(l):
            if (i > l or j >= k) and v.has_key((i,j,k,l)):
              lines.append("%20.16f%4d%4d%4d%4d" % (v[(i,j,k,l)], i+1,j+1,l+1,k+1))
              del v[(i,j,k,l)]
  assert(len(v) == 0)
  lines.append(empty_line())

def insert_symm_matrix(lines, v, nsites, seperation = True):
  assert(np.allclose(v.T - v, 0.))
  for i in range(nsites):
    for j in range(i+1):
      lines.append("%20.16f%4d%4d%4d%4d" % (v[i, j], i+1,j+1,0,0))
  if seperation:
    lines.append(empty_line())

def insert_matrix(lines, v, nsites):
  for i in range(nsites):
    for j in range(nsites):
      lines.append("%20.16f%4d%4d%4d%4d" % (v[i, j], i+1,j+1,0,0))
  lines.append(empty_line())

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

def get_block_result(orbitaltxt, nsites, nelec = -1, UHFB = False, p = False):
  if not os.path.exists(TmpDir):
    try:
      os.makedirs(TmpDir)
    except OSError:
      pass

  path = mkdtemp(prefix = "BLOCK", dir = TmpDir)
  if p:
    print "Temporary File Location"
    print "%s:%s" % (getoutput("echo $HOSTNAME"), path)
    print
  
  cwd = os.getcwd()
  os.chdir(path)
  WriteFile("DMETDUMP", orbitaltxt)

  # prepare "dmrg.conf" file
  configfile = []
  if nelec == -1:
    configfile.append("nelec %d\nspin 0\nhf_occ integral" % (nsites*2))
    configfile.append("bogoliubov" % (nsites*2))
  else:
    configfile.append("nelec %d\nspin 0\nhf_occ integral" % nelec)

  configfile.append(block_schedule(param_block["M"]))

  configfile.append("orbitals DMETDUMP")
  configfile.append("nonspinadapted\nonepdm")
  configfile.append("noreorder")
  configfile.append("outputlevel -1")

  configtxt = "\n".join(configfile) + "\n"

  WriteFile("dmrg.conf", configtxt)
  #print "Temporary Dir %s" % path
  cmd = ["mpiexec -np %d" % param_block["nproc"]]
  if param_block["bind"]:
    cmd.append("--bind-to-socket")
  cmd.append(param_block["exec"])
  cmd.append("dmrg.conf")
  cmd = " ".join(cmd)
  # call block
  blockoutput = getoutput(cmd)
  if p:
    print blockoutput
    print
  # read energy
  file_e = open("dmrg.e", "rb")
  energy = struct.unpack('d', file_e.read(8))[0]
  # read rdm and kappa
  if UHFB:
    rdm = readrdm("onepdm.0.0.txt", True)
  else:
    rdm = readrdm("spatial_onepdm.0.0.txt", False) / 2
  kappa = readrdm("spatial_pairmat.0.0.txt", False)
  if not UHFB:
    kappa = (kappa + kappa.T) / 2
  os.chdir(cwd)
  rmtree(path)
  return energy, rdm, kappa

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

def run_emb_bcs(h0, d0, vccdd, vcccd, vcccc, UHFB = False, p = False):
  if UHFB:
    nsites = h0[0].shape[0]
    orbitalfile = []
    orbitalfile.append(" &BCS NORB=%2d," % nsites)
    orbitalfile.append("  ORBSYM=%s" % ("1,"*nsites))
    orbitalfile.append("  ISYM=1,")
    orbitalfile.append("  IUHF=1,")
    orbitalfile.append(" &END")
    # vccdd_alpha, vccdd_beta, vccdd_ab
    for v in vccdd:
      insert_ccdd(orbitalfile, v, nsites)
    # vcccd_a, vcccd_b
    for v in vcccd:
      insert_cccd(orbitalfile, v, nsites)
    # vcccc
    insert_cccc(orbitalfile, vcccc, nsites, True)
    # h0a, h0b
    for h in h0:
      insert_symm_matrix(orbitalfile, h, nsites)
    # d0
    insert_matrix(orbitalfile, d0, nsites)
 
  else:
    nsites = h0.shape[0]
    orbitalfile = []
    orbitalfile.append(" &BCS NORB=%2d," % nsites)
    orbitalfile.append("  ORBSYM=%s" % ("1,"*nsites))
    orbitalfile.append("  ISYM=1,")
    orbitalfile.append(" &END")
    # vccdd
    insert_ccdd(orbitalfile, vccdd, nsites)
    # vcccd (a)
    insert_cccd(orbitalfile, vcccd, nsites)
    # vcccc
    insert_cccc(orbitalfile, vcccc, nsites, False)
    # h0
    insert_symm_matrix(orbitalfile, h0, nsites)
    # d0
    insert_symm_matrix(orbitalfile, d0, nsites)
  
  orbitalfile.append(empty_line()) # core energy is zero
  
  orbitaltxt = "\n".join(orbitalfile) + "\n"

  return get_block_result(orbitaltxt, nsites, UHFB, p)

def run_emb_normal(h0, vccdd, UHF = False, p = False):
  if UHF:
    nsites = h0[0].shape[0]
    orbitalfile = []
    orbitalfile.append(" &FCI NORB=%2d," % nsites)
    orbitalfile.append("  ORBSYM=%s" % ("1,"*nsites))
    orbitalfile.append("  ISYM=1,")
    orbitalfile.append("  IUHF=1,")
    orbitalfile.append(" &END")
    # vccdd_alpha, vccdd_beta, vccdd_ab
    for v in vccdd:
      insert_ccdd(orbitalfile, v, nsites)
    # h0a, h0b
    for h in h0:
      insert_symm_matrix(orbitalfile, h, nsites)
  else:
    nsites = h0.shape[0]
    orbitalfile = []
    orbitalfile.append(" &FCI NORB=%2d," % nsites)
    orbitalfile.append("  ORBSYM=%s" % ("1,"*nsites))
    orbitalfile.append("  ISYM=1,")
    orbitalfile.append(" &END")
    # vccdd
    insert_ccdd(orbitalfile, vccdd, nsites, false)
    # h0
    insert_symm_matrix(orbitalfile, h0, nsites, false)
  
  orbitalfile.append(empty_line()) # core energy is zero
  
  orbitaltxt = "\n".join(orbitalfile) + "\n"

  return get_block_result(orbitaltxt, nsites, UHF, p)
