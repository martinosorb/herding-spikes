import numpy as np
import pylab
import matplotlib as mplt
import h5py
import SpkD_v25

###Example script to use the spike detection in forward and reverse direction and merge the detected events
Folder='/home/xxx'
TxtFolderForward='/TxtFilesForward/'#where the .txt files from the forward detection are
TxtFolderReverse='/TxtFilesReverse/'#where the .txt files from the reverse detection are
HdfFolder='/HdfFiles/'#where the .hdf files should go
PlotFolder='/Plots/'#folder for plots,
Files=('chip82_rec00_19_10_2011_basal',)#filenames (without '_XXX.txt')

iii=0
PlotFile=Folder+Subfolder+PlotFolder+Files[iii]
TxtFileForward=Folder+Subfolder+TxtFolderForward+Files[iii]
TxtFileReverse=Folder+Subfolder+TxtFolderReverse+Files[iii]
HdfFileForward=Folder+Subfolder+HdfFolder+Files[iii]+'_forward_v25.hdf5'
HdfFileReverse=Folder+Subfolder+HdfFolder+Files[iii]+'_reverse_v25.hdf5'
HdfFileMerged=Folder+Subfolder+HdfFolder+Files[iii]+'_merged_v25.hdf5'
if removeCh==-1:
	NoisyChFile=Folder+Subfolder+HdfFolder+Files[iii]+'_NoisyChannels.hdf5'
else:
	NoisyChFile=''

NCh, tMax, Sampling=SpkD_v25.readInfoFile(TxtFileForward, HdfFile)
SpkD_v25.readAvgFile(TxtFileForward, HdfFileForward, NCh)
NSpk=SpkD_v25.readSpikesFile(TxtFileForward, HdfFileForward, NoisyChFile, NCh, removeCh, tMax)
SpkD_v25.readShapesFile(TxtFileForward, HdfFileForward, NSpk, Sampling, Lspike=4)
NCh, tMax, Sampling=SpkD_v25.readInfoFile(TxtFileReverse, HdfFileReverse)
SpkD_v25.readAvgFile(TxtFileReverse, HdfFileReverse, NCh)
NSpk=SpkD_v25.readSpikesFile(TxtFileReverse, HdfFileReverse, NoisyChFile, NCh, removeCh, tMax)
SpkD_v25.readShapesFile(TxtFileReverse, HdfFileReverse, NSpk, Sampling, Lspike=4)
SpkD_v25.MergeForwardBackward(HdfFileForward, HdfFileReverse, HdfFileMerged\
, IncludeLongSpikesForward=False, IncludeLongSpikesBackward=True, DFrames=3, MaxDist=1.)

