from inputs import Input, dump_inputfile
import basicForNormal
import block_iface
import helpers
import geometry
import normalDmet
import results
import typeclass
try:
    from main import *
except:
    def main(*args):
        print 'dmet is not running'
