#import pyximport
#pyximport.install()

from detect import detect
import h5py
import numpy as np

# raw data file
rawpath = 'data/'
rawfile = rawpath+'P29_16_05_14_retina02_left_stim2_smallarray_fullfield_raw3'
print(rawfile)

# get sampling rate
# path = '/data/MEA/LightStim/P29_16_07_14/'
# hpath = '/HdfFilesSpkD45_dev/'
# spikeFile = path + hpath + 'P29_16_05_14_retina02_left_stim2_smallarray_fullfield_v28_clustered_0.30_0.28_align'
# sf = h5py.File(spikeFile + '.hdf5' , 'r')
# sampling = sf['Sampling'].value
# print(sampling)
# sf.close()
sampling = 23199.0903585

# run detection
detect(rawfile, sampling)
