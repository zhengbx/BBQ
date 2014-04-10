import numpy as np
import numpy.linalg as la

class FFragment(object):
  def __init__(self, sites, emb_method, fitting_method):
    self.sites = np.array(sites)
    self.method = emb_method
    self.fitting = fitting_method

  def __str__(self):
    s = "Sites: %s   Impurity Solver: %s   Fitting Method: %s" % (self.sites, self.method, self.fitting)
    return s

  def get_sites(self):
    return self.sites

  def get_emb_method(self):
    return self.method
