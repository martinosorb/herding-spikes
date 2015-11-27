# Post-processing of detected spikes

Developed and written by [Oliver Muthmann](ollimuh@googlemail.com).

Reference: J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.


## Spikes detected with the online method

Methods used for the analysis are defined in SpkD_online.py. These are called by a script (Analysis_online.py), where parameters and filenames need to be set.

SpkD_online.readSpikesFile(TxtFile, HdfFile, NoisyChFile='', recCh=4096, removeCh=0, tMax=1200, Sampling=7563):

First, text files are read and converted into HDF5 format (under 'RawEvents/...'), along with parameters used during the spike detection (number of recorded channels, duration of the recording, sampling rate). (Optionally, a .hdf file containing the indices of noisy channels can be provided to ignore events in those channels.)

SpkD_online.IsolatedSpikes(HdfFile, DFrames=2, MaxDist=1.5)

To ensure that spikes are not detected on multiple electrodes, we remove events where an event of higher amplitude can be found within DFrames (in units of frames) temporal and MaxDist (in units of electrode distances) spatial distance. This creates arrays of time stamps ('Times'), event amplitudes ('Amplitudes') and channel ids ('Channels') of isolated events in the .hdf file. Additionally to channel ids (which are not used for the interpolating detection), event locations in (y,x) coordinates (electrodes centered at {0.5,...,63.5}) are saved as 'Locations'.

## Spikes detected with the interpolation method

Methods used for the analysis are defined in SpkD_interp.py. These are called by a script (Analysis_interp.py), where parameters and filenames need to be set.

SpkD_interp.readInfoFile(TxtFile, HdfFile):

Here parameters are already saved in the text files and do not need to be given manually. This creates a .hdf file and writes all the parameters of the detection.

SpkD_interp.readAvgFile(TxtFile, HdfFile):

This simply writes the global voltage fluctuations (as 'GlobalVoltageFluctuations/medianVoltage') into the .hdf file.

SpkD_interp.readSpikesFile(TxtFile, HdfFile, NoisyChFile, NCh, removeCh, tMax):

This writes the detected spikes, timestamps and channel ids to the .hdf file. Further, it returns the number of detected events to obtain the appropriate dimensions of arrays in the .hdf file.

SpkD_interp.readShapesFile(TxtFile, HdfFile, NSpk, Sampling, Lspike=5):

Localization of spikes and writing out averaged shapes. Parameters here are the number of events detected, the sampling rate, and the length of the temporal window in which the peak of the spike is expected.

SpkD_interp.IsolatedSpikes(HdfFile, IncludeLongSpikes=False, DFrames=2, MaxDist=1.):

Only retain events with largest amplitudes in a close spatiotemporal surrounding. This creates arrays of time stamps ('Times'), event amplitudes ('Amplitudes'), spatially averaged shapes ('Shapes') and locations ('Locations', in (y,x) coordinates) of isolated events in the .hdf file. The interpolating spike detection will detect events even if they do not repolarize. However, it appears that non repolarizing events contain a higher amount of false positives. For a more conservative detection, such 'long' events ('RepolarizingSpikes' is False) can be excluded from the analysis.

## Data plotting

A couple of elementary methods for visualizing the results are provided in SpkD_plot.py. These are

SpkD_plot.Rasterplot Creates colour coded rasters. A projection of the spatial locations on the x- or y-axis is plotted against time, and the colour can be used for any other property of events (location, amplitude, probability...).

SpkD_plot.Scatterplot Creates a colour coded spatial scatterplot of the quantity of interest.

SpkD_plot.Densityplot Creates a spatial histogram.

SpkD_plot.Matrixplot To visualize matrices.

SpkD_plot.Clusterplot Visualizes clusters from the clustering method.

SpkD_plot.PeriEventActivityPlot Plots a spatial map of biases in the activity of surrounding events.


## Correlation analysis

The correlation analysis requires longer datasets in order to detect correlations. Therefore raw data files for that would be unreasonably large in size and we provide reduced data in .hdf format where the detection and localization was already performed.

The first step determines a correlation index for each event. Parameters are half of the length of the window (in ranks) to detect correlations (NextN), the frequency (events per hour and bin) of added Poisson events (NPoissonNoise), the interval (in frames) between clock spikes (dClockspikes) and the number of largest local maxima to discard as reference channels as they might correspond to noisy channels (Nignore).

For the interpolating detection, additional parameters are the spatial resolution (Res, in bins per 42µm), the spatial resolution for determination of local maxima (RefRes) and the number of frames to discard around recalibration events when determining correlations (removeRecalib). Further, there is a possibility to adjust the significance level for correlations (probC). This step writes correlation indices for each event ('Corr') to the .hdf file.

SpkD_online.CorrelationAnalysis1(HdfFile, NextN=100, NPoissonNoise=360,\ dClockspikes=4, Nignore=10)

SpkD_interp.CorrelationAnalysis1(HdfFile, NextN=100, NPoissonNoise=40,\ dClockspikes=4, removeRecalib=3, Nignore=10, probC=0.01, Res=3, RefRes=5)

The second step performs a comparison between surrogate and detected events and determines a fraction of false positives (saved in the .hdf file as 'CorrelationAnalysis/Noise') and a probability estimate for each event to correspond to a real spike (saved as 'CorrelationAnalysis/Probability'). It also stores the p-values of the KS-test ('CorrelationAnalysis/pValue')

SpkD_online.CorrelationAnalysis2(HdfFile)

SpkD_interp.CorrelationAnalysis2(HdfFile)
