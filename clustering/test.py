from herdingspikes import spikeclass, QualityMeasures
import numpy as np
import numpy.testing as npt
import unittest


class TestHS(unittest.TestCase):
    n1 = 10000
    data1 = np.random.normal(size=(2, n1))
    data2 = np.random.normal(3, 1, size=(2, n1))
    data = np.hstack([data1, data2])
    n = 2*n1
    princcomp = np.random.random((2, n))

    def test_ms_simple(self):
        # the simplest import
        O = spikeclass(self.data)
        self.assertEqual(O.NData(), self.n)
        # combined mean shift, expected result 2 clusters
        O.CombinedMeanShift(2., 0.001, PrincComp=self.princcomp, njobs=1)
        self.assertEqual(O.NClusters(), 2)
        # check ClusterSizes and the cluster sizes
        npt.assert_array_almost_equal(O.ClusterSizes()/self.n,
                                      [0.5, 0.5], decimal=3)

    def test_gaussian_overlap(self):
        O = spikeclass(self.data)
        O.CombinedMeanShift(2., 0.001, PrincComp=self.princcomp, njobs=1)
        fakePCA = np.empty_like(self.data)
        QM = QualityMeasures(O, scorePCA=fakePCA)
        data = []
        for ind in [0, 1]:
            data.append(O.Locations()[:, O.ClusterID() == ind].T)
        estCent = np.array([np.mean(d, axis=0) for d in data])
        g1 = QM._fit_gaussian_mixture(data)
        npt.assert_array_almost_equal(g1.means_, estCent, decimal=1)
        g2 = QM._fit_gaussian_individuals(data)
        npt.assert_array_almost_equal(g2.means_, estCent, decimal=1)
        npt.assert_array_almost_equal(g1.covars_, g2.covars_, decimal=1)
        npt.assert_array_almost_equal(g1.weights_, g2.weights_, decimal=2)

        m1 = QM._data_gaussian_overlap(data, fit_mode="mixture")
        m2 = QM._data_gaussian_overlap(data, fit_mode="individuals")
        npt.assert_array_almost_equal(m1, m2, decimal=2)


if __name__ == '__main__':
    unittest.main()
