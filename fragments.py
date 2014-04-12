#
# File: fragments.py
# Author: Bo-Xiao Zheng <boxiao@princeton.edu>
#

import numpy as np
import numpy.linalg as la

class FFragment(object):
    def __init__(self, sites, emb_method, fitting_method, name = None):
        self.sites = np.array(sites)
        self.method = emb_method
        self.fitting = fitting_method
        self.name = name

    def __str__(self):
        s = "Sites: %20s      Impurity Solver: %10s       Fitting Method: %10s" % (self.sites, self.method, self.fitting)
        if self.name is not None:
            s = "Name:%8s " % self.name + s
        return s

    def get_sites(self):
        return self.sites

    def get_emb_method(self):
        return self.method
