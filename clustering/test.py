from herdingspikes import spikeclass
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

if __name__ == '__main__':
    unittest.main()
