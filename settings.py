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

#FciExecutable = "/home/cgk/dev/fci/fci"
FciExecutable = "/home/qiaoni/Documents/work/dev/dmet-hubbard/fci.20121221/fci"
param_block = {
  "exec": "/home/boxiao/block/block.spin_adapted",
  "M": 400,
  "nproc": 16,
  "bind": False,
}
CCExecutable = None
#FciExecutable = "/home/qiaoni/Documents/work/dev/dmet-hubbard/fci.20121221/fci"

TmpDir = None
TmpDir = "./fci"
LibDir = "./libdir"
#TmpDir = "/home/qiaoni/Templates/"

#TmpDir = "/dev/shm/"
# ^- if not None, use this as prefix for temporary files. Otherwise
#    use whatever mkdtmp() generates (likely related to TMPDIR environment
#    variable)
