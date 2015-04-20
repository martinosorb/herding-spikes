import numpy as np
import pylab
import os.path
import matplotlib as mplt
import h5py
import SpkD_interp
import SpkD_plot

###This is the program that writes the .txt Files into a .hdf Files
# and estimates the spatial origins of spikes, marks events detected in multiple channels,

#...where to find/save the files
TxtFile='BW_Chip_167_PT_Prep_12_11_30_Exp_24_12_2012_2500msec_burstActivity'
HdfFile=TxtFile+'_SpkD_interp.hdf5'
PlotFile='Chip167_Plots'

#--------------------------------------------------------------------
###Parameters
#which channels to remove
removeCh=0#-1: using a list from a file, 0: no removal, >0 cutoff frequency in Hz
NoisyChFile=''
#Localization
NCh, tMax, Sampling=SpkD_interp.readInfoFile(TxtFile, HdfFile)
SpkD_interp.readAvgFile(TxtFile, HdfFile, NCh)
NSpk=SpkD_interp.readSpikesFile(TxtFile, HdfFile, NoisyChFile, NCh, removeCh, tMax)
SpkD_interp.readShapesFile(TxtFile, HdfFile, NSpk, Sampling, Lspike=5)
SpkD_interp.IsolatedSpikes(HdfFile, IncludeLongSpikes=False, DFrames=3, MaxDist=1.)

#Plotting
###Rasters
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'_Raster.png', FieldName=''\
, tMin=0, tMax=2.5, yMin=0, yMax=64, cMin=0, cMax=64, cRes=20\
, RepolSpikes=True, LongSpikes=True, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.2)
#Amplitudes
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'_Raster_amplitudes.png'\
, FieldName='Amplitudes'\
, tMin=0, tMax=2.5, yMin=0, yMax=64, cMin=0, cMax=10., cRes=20\
, RepolSpikes=True, LongSpikes=False, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.2)
##Scatterplot all spikes amplitudes
SpkD_plot.Scatterplot(HdfFile, FileName=PlotFile+'_Scatter_repolarizing.png', FieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=10, cRes=20\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.1)
###Densityplot
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'_Density_repolarizing.png', FieldName=''\
, ProbFieldName='', xMin=0, xMax=64, yMin=0, yMax=64, cMin=10, cMax=10000, Res=8, logScale=True\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
#Amplitudes
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'_Density_amplitudes.png', FieldName='Amplitudes'\
, ProbFieldName='', xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=10, Res=8, logScale=False\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
###Global voltage fluctuations
SpkD_plot.GlobalVoltageplot(HdfFile, FileName=PlotFile+'_GlobalVCorr.png'\
, cMin=-0.2, cMax=0.4, CorrCmap=pylab.cm.hsv)
