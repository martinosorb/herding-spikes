import unittest
from spikeclus import spikeclass,ReadInterpolated
import matplotlib.pyplot as plt

class ClusterTest(unittest.TestCase):
    def test_import_data(self):
        O = spikeclass('saved-fortest.hdf5')
        self.assertEqual(O.NData(),160)
    def test_clustering_fsc(self):
        O = spikeclass('saved-fortest.hdf5')
        O.skMS()
        self.assertEqual(O.NClusters(),2)
        a = O.FilterSmallClusters(20)
        self.assertEqual(O.NClusters(),1)
        self.assertTrue(all(a == range(150)))
        self.assertTrue(all(O.ClusterID() == 0))
    def test_fld(self):
        O = spikeclass('saved-fortest.hdf5')
        O.FilterLowDensity(1,[10,10])
        self.assertEqual(O.NData(),140)
    def test_shapes(self):
        O = ReadInterpolated('Retina02_LeftEye_Stim01_whitenoise100ms_SpkD45_v18.hdf5')
        O.Crop(1.5,4.5,56,59)
        O.ShapePCA(2)
        
if __name__ == '__main__':
    unittest.main()