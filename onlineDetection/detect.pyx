# distutils: language = c++
# distutils: sources = SpkDonline.cpp

import cython
import numpy as np
cimport numpy as np
cimport cython
import h5py
from ctypes import CDLL
import ctypes


cdef extern from "SpkDonline.h" namespace "SpkDonline":
    cdef cppclass Detection:
        Detection() except +
        void InitDetection(long nFrames, double nSec, int sf, int NCh, long ti, long int * Indices)
        void SetInitialParams(int thres, int maa, int ahpthr, int maxsl, int minsl)
        void openSpikeFile(const char * name)
        void MedianVoltage(unsigned short * vm)
        void MeanVoltage(unsigned short * vm)
        void Iterate(unsigned short * vm, long t0)
        void FinishDetection()


def detect(rawfilename, sfd):
    rf = h5py.File(rawfilename + '.hdf5', 'r')
    cx, cy = rf['x'].value, rf['y'].value
    x1, x2, y1, y2 = np.min(cx), np.max(cx), np.min(cy), np.max(cy)
    nRecCh = len(cx)
    nFrames = len(rf['Ch0'])
    nSec = nFrames / sfd  # the duration in seconds of the recording
    sf = int(sfd)
    nDumpFrames = sf * 5  # nFrames; # how many frames to analyze
    nSec = nDumpFrames / sfd

    tInc = 10000
    tCut = int(1.0 * sf / 1000 + 1.0 * sf / 1000 + 6)

    print("# Sampling rate: " + str(sf))
    print("# Number of recorded channels: " + str(nRecCh))
    print("# Analysing frames: " + str(nDumpFrames) + ", Seconds:" + str(nSec) + ", tCut:" + str(tCut) + ", tInc:" + str(tInc))

    # Messy! To be consistent, X and Y have to be swappped
    cdef np.ndarray[long, mode = "c"] Indices = np.zeros(nRecCh, dtype=long)
    for i in range(nRecCh):
        Indices[i] = (cy[i] - 1) + 64 * (cx[i] - 1)

    cdef Detection * det = new Detection()
    det.InitDetection(nFrames, nSec, sf, nRecCh, tInc, &Indices[0])
    det.SetInitialParams(9, 5, 0, 8, 3)

    # open output file
    spikefilename = str.encode(rawfilename + "_Spikes.txt")
    det.openSpikeFile(spikefilename)

    cdef np.ndarray[unsigned short, ndim = 1, mode="c"] vm = np.zeros((nRecCh*tInc), dtype=ctypes.c_ushort)
    for t0 in range(0, nDumpFrames - tInc, tInc - tCut):
        for c in range(nRecCh):
            vm[c*tInc:c*tInc+tInc] = rf['Ch'+str(c)][t0:t0+tInc].astype(dtype=ctypes.c_ushort)
        #det.MedianVoltage(&vm[0]) # only use if there are many channels
        det.MeanVoltage(&vm[0]) # use if there are few channels
        det.Iterate(&vm[0], t0)
    det.FinishDetection ();
