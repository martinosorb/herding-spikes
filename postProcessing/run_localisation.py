import os.path
import sys
import time
import pylab
import SpkD_v28
import SpkD_plot_v24

# Parameters
# which channels to remove
# -1: using a list from a file, 0: no removal, >0 cutoff
removeCh = 0

Folder, Filename = os.path.split(sys.argv[1])

TxtFile = Folder + '/' + Filename.replace('.brw', '_INT')
print("Input: " + TxtFile)
HdfFile = Folder + '/' + Filename.replace('.brw', '') + '_v28.hdf5'
print("Output: " + HdfFile)
if removeCh == -1:
    NoisyChFile = Folder + '/' + \
        Filename.replace('.brw', '') + '_NoisyChannels.hdf5'
    print("Noisy channels are in: " + NoisyChFile)
else:
    NoisyChFile = ''

start_time = time.time()
NCh, tMax, Sampling = SpkD_v28.readInfoFile(TxtFile, HdfFile)
SpkD_v28.readAvgFile(TxtFile, HdfFile, NCh)
NSpk = SpkD_v28.readSpikesFile(
    TxtFile, HdfFile, NoisyChFile, NCh, removeCh, tMax)
SpkD_v28.readShapesFile(TxtFile, HdfFile, NSpk)
SpkD_v28.IsolatedSpikes(
    HdfFile, IncludeLongSpikes=False, DFrames=3, MaxDist=1.)
print("time taken: %s sec" % (time.time() - start_time))

print("generating plots")
PlotFile = TxtFile
# Rasterplot
SpkD_plot_v24.Rasterplot(HdfFile,
                         FileName=PlotFile + '_Raster_repolarizing.png',
                         FieldName='', tMin=0, tMax=20, yMin=0, yMax=64,
                         cMin=0, cMax=64, cRes=20, RepolSpikes=True,
                         LongSpikes=False, axis=0, ampThreshold=0,
                         pThreshold=0, bgColor=pylab.cm.gray(0.5),
                         SpikesCmap=pylab.cm.hot_r,
                         LongSpikesCmap=pylab.cm.cool, alpha=0.2)
# Amplitudes
SpkD_plot_v24.Rasterplot(HdfFile,
                         FileName=PlotFile + '_Raster_amplitudes.png',
                         FieldName='Amplitudes', tMin=0, tMax=20, yMin=0,
                         yMax=64, cMin=0, cMax=10., cRes=20,
                         RepolSpikes=True, LongSpikes=False, axis=0,
                         ampThreshold=0, pThreshold=0,
                         bgColor=pylab.cm.gray(0.5),
                         SpikesCmap=pylab.cm.hot_r,
                         LongSpikesCmap=pylab.cm.cool, alpha=0.2)

# Scatterplot
SpkD_plot_v24.Scatterplot(HdfFile,
                          FileName=PlotFile + '_Scatter_repolarizing.png',
                          FieldName='Amplitudes', xMin=0, xMax=64, yMin=0,
                          yMax=64, cMin=0, cMax=10, cRes=20,
                          RepolSpikes=True, LongSpikes=False,
                          ampThreshold=0, pThreshold=0,
                          bgColor=pylab.cm.gray(0.5),
                          SpikesCmap=pylab.cm.hot_r,
                          LongSpikesCmap=pylab.cm.cool, alpha=0.1)

# Densityplot
SpkD_plot_v24.Densityplot(HdfFile,
                          FileName=PlotFile + '_Density_repolarizing.png',
                          FieldName='', xMin=0, xMax=64, yMin=0, yMax=64,
                          cMin=10, cMax=10000, Res=8, logScale=True,
                          RepolSpikes=True, LongSpikes=False,
                          ampThreshold=0, pThreshold=0,
                          bgColor=pylab.cm.gray(0.5),
                          SpikesCmap=pylab.cm.hot_r, alpha=1.)

# Amplitudes
SpkD_plot_v24.Densityplot(HdfFile,
                          FileName=PlotFile + '_Density_amplitudes.png',
                          FieldName='Amplitudes', xMin=0, xMax=64, yMin=0,
                          yMax=64, cMin=0, cMax=10, Res=8, logScale=False,
                          RepolSpikes=True, LongSpikes=False,
                          ampThreshold=0, pThreshold=0,
                          bgColor=pylab.cm.gray(0.5),
                          SpikesCmap=pylab.cm.hot_r, alpha=1.)
