import os.path
import SpkD_v28

# This is an example script that writes the .txt Files into a .hdf Files
# and estimates the spatial origins of spikes and marks events detected in
# multiple channels.

# ...where I had the files on my computer
Folder = '/home/muthmann/backup/SpikeAnalysis/'
Subfolder = 'SpkD_v26/'
TxtFolder = '/TxtFiles/'  # where the .txt files are
HdfFolder = '/HdfFiles/'  # where the .hdf files should go
Files = ('P10_6Feb15_ret2_waves_raw',)  # filenames (without '_XXX.txt')
# --------------------------------------------------------------------
# Parameters
# which channels to remove
removeCh = 0  # -1: using a list from a file, 0: no removal, >0 cutoff
#               frequency in Hz empty chip recording
# --------------------------------------------------------------------------
if not os.path.exists(Folder+Subfolder+HdfFolder):
    os.mkdir(Folder+Subfolder+HdfFolder)
for iii in range(len(Files)):
    print(Files[iii])
    TxtFile = Folder+Subfolder+TxtFolder+Files[iii]
    HdfFile = Folder+Subfolder+HdfFolder+Files[iii]+'_v28.hdf5'
    if removeCh == -1:
        NoisyChFile = Folder+Subfolder+HdfFolder+Files[iii] + \
                        '_NoisyChannels.hdf5'
    else:
        NoisyChFile = ''
    NCh, tMax, Sampling = SpkD_v28.readInfoFile(TxtFile, HdfFile)
    SpkD_v28.readAvgFile(TxtFile, HdfFile, NCh)
    NSpk = SpkD_v28.readSpikesFile(TxtFile, HdfFile, NoisyChFile, NCh,
                                   removeCh, tMax)
    SpkD_v28.readShapesFile(TxtFile, HdfFile, NSpk)
    SpkD_v28.IsolatedSpikes(HdfFile, IncludeLongSpikes=False,
                            DFrames=3, MaxDist=1.)
