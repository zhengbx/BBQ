#!/usr/bin/env python

import unittest
import numpy
import BBQ


def rand_hermit_mat(n, seed=1):
    numpy.random.seed(seed)
    a = numpy.random.random((n,n))
    return a + a.T

class KnowValues(unittest.TestCase):
    def test_resmin(self):
        print 'todo: test resmin'

    def test_FitCorrelationPotential(self):
        print 'todo: test FitCorrelationPotential'

    def test_MakeEmbeddingBasis(self):
        norb = 8
        nocc = 3
        fock = rand_hermit_mat(norb, 15)
        e, c = numpy.linalg.eigh(fock)
        dm = numpy.dot(c[:,:nocc],c[:,:nocc].T)
        impsites = [0, 1]
        b1 = BBQ.basicForNormal.MakeEmbeddingBasis(impsites, dm, 1e-8)
        b2 = BBQ.basicForNormal.MakeEmbeddingBasis1(impsites, dm, 1e-8)
        self.assertAlmostEqual(abs(b1-b2).sum(), 0., 14)

    def test_fci_solver_r(self):
        m = BBQ.normalDmet.NormalDmet('RHF',nElec=6)
        fock = rand_hermit_mat(6, 15)
        EmbRdm = numpy.eye(6)
        res = m.ImpSolver(fock, EmbRdm, 'Fci', U=1)
        self.assertAlmostEqual(res.Energy, -1.0712572404440699, 10)
        d = numpy.array((0.004944545542672,
                         0.013674600980115,
                         0.034250514010321,
                         1.974467334084587,
                         1.978031303849934,
                         1.994631701532369,))
        diff = abs(numpy.linalg.eigh(res.Rdm)[0] - d).sum()
        self.assertAlmostEqual(diff, 0., 10)

    def test_fci_solver_u(self):
        m = BBQ.normalDmet.NormalDmet('UHF', nElecA=3, nElecB=4)
        fock = numpy.zeros((12,12))
        fock[ ::2, ::2] = rand_hermit_mat(6, 14)
        fock[1::2,1::2] = rand_hermit_mat(6, 15)
        EmbRdm = numpy.eye(6)
        EmbRdm[0,0] += 1
        res = m.ImpSolver(fock, EmbRdm, 'Fci', U=1)
        self.assertAlmostEqual(res.Energy, 17.02787257111538, 10)
        d = numpy.array((0.001597212773413,
                         0.004283739517800,
                         0.004620116081424,
                         0.005536223219325,
                         0.006816761844207,
                         0.992183999123251,
                         0.993298137788704,
                         0.996709467732277,
                         0.996811438769326,
                         0.998783976376671,
                         0.999587386596884,
                         0.999771540176729))
        diff = abs(numpy.linalg.eigh(res.Rdm)[0] - d).sum()
        self.assertAlmostEqual(diff, 0., 10)

    def test_main(self):
        print 'todo: test main'


if __name__ == "__main__":
    print "Full Tests for Lattice DMET"
    unittest.main()
