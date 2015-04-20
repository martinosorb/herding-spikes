import numpy as np
import pylab
import os.path
import matplotlib as mplt
import h5py
import SpkD_online
import SpkD_plot

###This is the program that writes the .txt File of the online version into a .hdf File,
# marks events detected in multiple channels,
# and runs a correlation analysis.

#...where to find/save the files
TxtFile='chip82_rec00_19_10_2011_basal'
HdfFile=TxtFile+'_SpkD_online.hdf5'
PlotFile='Chip82_online_Plots'

#--------------------------------------------------------------------
###Parameters
#which channels to remove
removeCh=0#-1: using a list from a file, 0: no removal, >0 cutoff frequency in Hz
NoisyChFile=''

#BW_Chip_192_lowDense_PT_phase_00
#empty chip recording
#P11_14May12_spont_day3_9
#chip82_rec00_19_10_2011_basal
#
#--------------------------------------------------------------------------
#convert to .hdf
NSpk=SpkD_online.readSpikesFile(TxtFile, HdfFile, NoisyChFile, recCh=4096, removeCh=0, tMax=1200, Sampling=7563)
SpkD_online.IsolatedSpikes(HdfFile, DFrames=2, MaxDist=1.5)
#correlation analysis
SpkD_online.CorrelationAnalysis1(HdfFile, NextN=100, NPoissonNoise=360,\
dClockspikes=4, Nignore=10)
SpkD_online.CorrelationAnalysis2(HdfFile)

#Plotting
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'Raster.png', FieldName=''\
, tMin=0, tMax=20, yMin=0, yMax=64, cMin=0, cMax=64, cRes=20\
, RepolSpikes=True, LongSpikes=False, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.2, onlineVersion=True)
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'Raster_probabilities.png', FieldName='CorrelationAnalysis/Probability'\
, tMin=0, tMax=20, yMin=0, yMax=64, cMin=0, cMax=1., cRes=20, cLabel='probability'\
, RepolSpikes=True, LongSpikes=False, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.3, onlineVersion=True)
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'Density.png', FieldName=''\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=10, cMax=10000, Res=2, logScale=True\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'Density_probabilities.png', FieldName='CorrelationAnalysis/Probability'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., Res=2, logScale=True, cLabel='probability'\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
