import numpy as np
import pylab
import os.path
import matplotlib as mplt
import h5py
import SpkD_v28
import SpkD_plot_v24

###This is the program that
# runs a correlation analysis and, based on these results, does a clustering.


#...where to find/save the files
#HdfFile='chip82_rec00_19_10_2011_basal_SpkD_interp.hdf5'
HdfFile='P11_14May12_spont_day3_9_SpkD_interp.hdf5'
#HdfFile='empty chip recording_SpkD_interp.hdf5'
#PlotFile='Chip82_interp_Plots'
PlotFile='Retina_interp_Plots'
#NoisyChFile='chip82_rec00_19_10_2011_basal_NoisyChannels.hdf5'
NoisyChFile='P11_14May12_spont_day3_9_NoisyChannels.hdf5'
#NoisyChFile='empty chip recording_NoisyChannels.hdf5'
#--------------------------------------------------------------------

SpkD_interp.CorrelationAnalysis1(HdfFile, NextN=100, NPoissonNoise=40,\
dClockspikes=4, removeRecalib=3, Nignore=10, probC=0.01, Res=3, RefRes=5)
SpkD_interp.CorrelationAnalysis2(HdfFile)
SpkD_interp.Clustering(HdfFile, fNextMin=0.02,  fNext=0.1, Resolution=12,\
gradientThr=0.5, MaxDist=1., pOffset=0.2)
SpkD_interp.FindNoisyChannels([HdfFile,], NoisyChFile, nThreshold=3)
SpkD_interp.PeriEventActivity(HdfFile, NoisyChFile, minAmp=1.7, Res=12,Ns0=5./12,Na0=4, nNext=200)


#Plotting
###Raster
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'_Raster_repolarizing.png', FieldName=''\
, tMin=0, tMax=20, yMin=0, yMax=64, cMin=0, cMax=64, cRes=20\
, RepolSpikes=True, LongSpikes=False, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.2)
#Probabilities
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'_Raster_probabilities.png'\
, FieldName='CorrelationAnalysis/Probability'\
, tMin=0, tMax=20, yMin=0, yMax=64, cMin=0, cMax=1., cRes=20, cLabel='probability'\
, RepolSpikes=True, LongSpikes=False, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.2)
#Amplitudes
SpkD_plot.Rasterplot(HdfFile, FileName=PlotFile+'_Raster_amplitudes.png'\
, FieldName='Amplitudes'\
, tMin=0, tMax=20, yMin=0, yMax=64, cMin=0, cMax=10., cRes=20\
, RepolSpikes=True, LongSpikes=False, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.2)
###Scatterplot amplitudes
SpkD_plot.Scatterplot(HdfFile, FileName=PlotFile+'_Scatter_repolarizing.png', FieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=10, cRes=20\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.1)
#Probabilities
SpkD_plot.Scatterplot(HdfFile, FileName=PlotFile+'_Scatter_probabilities.png'\
, FieldName='CorrelationAnalysis/Probability'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., cRes=20, cLabel='probability'\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.1)
###Densityplot
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'_Density_repolarizing.png', FieldName=''\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=10, cMax=10000, Res=8, logScale=True\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
#Correlation indices
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'_Density_corrIndex.png'\
, FieldName='Corr'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., Res=8, logScale=False, cLabel='correlation index'\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
#Probabilities
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'_Density_probabilities.png'\
, FieldName='CorrelationAnalysis/Probability'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., Res=8, logScale=False, cLabel='probability'\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
#Amplitudes
SpkD_plot.Densityplot(HdfFile, FileName=PlotFile+'_Density_amplitudes.png', FieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=10, Res=8, logScale=False\
, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, alpha=1.)
#fraction of noise
SpkD_plot.Matrixplot(HdfFile, FileName=PlotFile+'_FracNoise.png', FieldName='CorrelationAnalysis/Noise'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., Res=3, logScale=False, cLabel='fraction of noise'\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r)
#p-Values
SpkD_plot.Matrixplot(HdfFile, FileName=PlotFile+'_pValue.png', FieldName='CorrelationAnalysis/pValue'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0.001, cMax=1., Res=3, logScale=True, cLabel='p-Value'\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r)
#Patterns
SpkD_plot.PeriEventActivityPlot(HdfFile, FileName=PlotFile+'_PeriEventActivity')
###Clusterplot
#SpkD_plot.Clusterplot(HdfFile, FileName=PlotFile+'_Cluster.png'\
#, xMin=0, xMax=64, yMin=0, yMax=64, cMin=10, cMax=10000, Res=12\
#, RepolSpikes=True, LongSpikes=False, ampThreshold=0, pThreshold=0, bgColor='w'\
#, ClusterCmap=pylab.cm.hsv)
