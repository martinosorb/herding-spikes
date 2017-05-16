import os.path
import sys
import time
from herdingspikes import spikeclass, ImportInterpolated

# Parameters for data clustering
h = 0.3  # kernel size
alpha = 0.28  # scaling of PCA relative to location

if __name__ == "__main__":
    Folder, Filename = os.path.split(os.path.abspath(sys.argv[1]))

    O = ImportInterpolated(sys.argv[1])
    scorePCA = O.ShapePCA(ncomp=2, white=True)
    O.CombinedMeanShift(h, alpha, PrincComp=scorePCA, mbf=10)
    print('Found ' + str(O.NClusters()) +
          ' clusters for ' + str(O.NData()) + ' spikes.')
    O.Save(Folder + '/' + Filename.replace('.hdf5',
                                           '_clustered_' + str(h) + '_' + str(alpha) + '.hdf5'))
