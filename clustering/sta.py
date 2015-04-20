"""
Copied by Martino from an IPython notebook on Wed 1 Apr 2015

@author: Sahar Pirmoradian
"""
import h5py
import numpy as np
import scipy.io
import os.path
from scipy import stats
import matplotlib as mplt
import matplotlib.pyplot as plt
import pylab
import scipy.optimize as opt

masterdr = 'non-git/sta'

# data file is an hdf5 file (spikes, channels, ...)
datafullfn = masterdr + '/hdf5/P38_06_03_14_retina02_lightstim_sorted.hdf5'

# stimulus file is a mat file that contains SamplingFrequency, frameSamples, frameTimes, pixels_allframes (see read_stimfile())  
stimfullfn = masterdr + '/mat/Retina02_LeftEye_Stim01_whitenoise100ms_stim.mat'

figpath = masterdr + '/results/' # where STA figures are stored
staparams_matfile = figpath + 'staparams' # it will be stored as a mat file
spksrc = 'su'  # spike source of our data, 'ch' (channels), 'cl' (clusters), 'su' (sorted units)
coord_spksrc = 'rc' # coordination of channles on the chip, 'rc' (row-column) or 'cr' (column-row) 
minamp = 7 # minimum amplitude of spikes
minspk = 30 # minimum spikes a cell must have in Full-field experiment to be considered in analyses; if the number of FF trials are 30, we assume a responding cell fires at least once to each trial, that is minimum spikes are 30
dur_sta = .5 # in sec, duration of STA 
dt_sta = 0.01 # in sec, temporal resolution of STA
dt_sta_graphics = .1 # in sec, graphically we show this temporal resolution on STA plots
spklatency_sta = .0  # always kept 0
figformat = '.png' # format of figures
max_var_acceptable = 20 # a fitted receptive field is acceptable if its center is fit with a variance that does not exceed this number
max_pval_acceptable = 1e-06 # p-value of STA peaks must not exceed this number
checkSTAPeakZFLG = True #if True, the STA peak is checked for significance
staplotFLG = True  #if True, STA plots are saved
PrintGaussParamFLG = False # if True, you can see the parameters of Gaussian fit
slctdchs_sta = np.array([13])# 'all' # if you like to compute STA of all cells, set slctdchs_sta = 'all'; otherwise if you compute only the STA of selected cells, e.g. 13 or 752, set slctdchs_sta = np.array([13,752])
plotJustRespUnitsFLG = True # if True, you only save the STA plots of cells with a significant STA peak
staSaveParamsFLG = True
no_gausfit_params = 7
chip_size = 2.67 # in mm
no_neighbpix = 729 # image size (27*27)
tableau = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229), 
              (0,0,128)] # color
#converting to rgb
tableau = [(c[0]/255., c[1]/255., c[2]/255.) for c in tableau]
    
# TODO (not used in this version)
figpath_rate = '.' # where Firing Rate figures are stored
figpath_mosaic = '.' # where figures of spatial arrangements of cells are stored
calcBiasidx4staFLG = False 
mosaicplotFLG=False
args_fr = None
mosaic_imgsize=729
plotJustFFRespUnitsFLG=False

def read_hdf5(fullfn, spksrc='ch'):
  # fullfn - string; full file name of an hdf5 file to be read
  # spksrc - 'ch' (channels), 'cl' (clusters), or 'su' (sorted units); spike source
  try:  
    myfile = h5py.File(fullfn,'r')
  except:
    raise IOError('\nI cannot open your hdf5 file: %s' % fullfn)

  spikes = np.array(myfile['Times'].value, dtype=int)
  amps = myfile['Amplitudes'].value
  freq = myfile['Sampling'].value
  loc2d, locspks2d = (np.array([]), np.array([]))
  if (spksrc == 'ch'):
    ch = myfile['Channels'].value
    loc = np.array([])
    nbins = 64
  elif (spksrc == 'cl'):
    locspks2d = myfile['Locations'].value
    ch = myfile['Cluster']['ClusterId'].value
    loc2d = myfile['Cluster']['CLocations'].value    
    loc = myfile['Cluster']['CLocationsS'].value  
    nbins = myfile['Cluster']['nBins'].value #768 #
    # removing clusters with id equal to -1
#    pdb.set_trace()
    valididxs = np.where(ch>=0)
    spikes = spikes[valididxs]
    amps = amps[valididxs]
    ch = ch[valididxs]    
    locspks2d = locspks2d[valididxs]    
       
  elif (spksrc == 'su'):
    ch = myfile['SortedUnits'].value
    markedch = myfile['MarkedChannels'].value
    valididxs = np.where(markedch>-1)
    spikes = spikes[valididxs]
    amps = amps[valididxs]
    ch = ch[valididxs]    
    loc = myfile['LocationsS'].value
    nbins = 64  
  data = {'spikes':spikes, 'amps': amps, 'ch':ch, 'freq':freq, 'nbins': nbins, 'loc':loc, 'loc2d':loc2d, 'locspks2d': locspks2d}
  return data


data = read_hdf5(datafullfn, spksrc=spksrc)
data

def read_stimfile(fullfn):
  # fullfn - string; a mat file containing the stimulus file whose variables are:
  # SamplingFrequency - float; sampling frequency of data acquisition, e.g. 7022
  # frameTimes - horizontal vector; times at which stimuli frames were presented, e.g. in a 15min white noise presentation with 100msec duration, 9000 images are presented, thus frameTimes is a <1x9000 double> vector  
  # frameSamples - horizontal vector; instead of time, it contains frame numbers at which stimuli were presented, i.e. frameSamples = frameTimes * SamplingFrequency 
  # pixels_allframes - 2D array; each row represents the pixels of the image presented at each frame, e.g. in a 15min white noise presentation with 100msec duration, 9000 images are presented, and if each image has 27*27 (=729) pixels, then pixels_allframes is a <9000x729 double> array  
  myfile = scipy.io.loadmat(fullfn)
  stimframes = np.array(myfile.get('frameSamples')[0], dtype='int')
  stimtimes = np.array(myfile.get('frameTimes')[0])
  stimpixels = np.array(myfile.get('pixels_allframes'), dtype='int')
  stimfreq = np.array(myfile.get('SamplingFrequency')[0][0])
  stim = {'stimframes':stimframes, 'stimtimes': stimtimes, 'stimfreq':stimfreq, 'stimpixels':stimpixels}  
  return stim

stim = read_stimfile(stimfullfn) 
stim

def get_spkidx(data, minamp=7, minframe=0, maxframe=1000):
  amps = data['amps']
  spikes = data['spikes']
  idxframe_min = np.where(spikes >= minframe)[0]
  elmidx_beforemaxf = next(x[0] for x in enumerate(idxframe_min) if spikes[x[1]] > maxframe)
  idxframe = idxframe_min[:elmidx_beforemaxf]
  idxamp = np.where(amps >= minamp)[0]
  idx = np.intersect1d(idxframe, idxamp)
  return idx

def get_hist_spkcnt(data, spkidx):
  ch = data['ch']
  no_ch=max(ch[spkidx]) + 1
  hist_spkcnt = np.histogram(ch[spkidx],bins=np.arange(no_ch+1))[0]
  return hist_spkcnt

def get_chlocation(ich, spksrc, nbins, loc, coord_spksrc='rc'):
  if 'array' in type(loc[ich]).__name__:
    locr = loc[ich][0]
    locc = loc[ich][1]
  elif spksrc == 'ch':
    locr = ich / nbins
    locc = np.mod(ich, nbins)
  elif spksrc == 'cl' or spksrc == 'su':
    if nbins:
      locr = loc[ich] / nbins
      locc = np.mod(loc[ich], nbins)
      locr = locr * 64 / nbins
      locc = locc * 64 / nbins  
    else:
      locr = loc[ich][0]
      locc = loc[ich][1]
  if coord_spksrc == 'cr':
    tmp = locr
    locr = locc
    locc = tmp
  return locr, locc

def get_neighborpixels(locr, locc, image_xsize, no_neighbpix_x):
  minr = locr - int(no_neighbpix_x/2)
  maxr = locr + int(no_neighbpix_x/2) + 1
  minc = locc - int(no_neighbpix_x/2)
  maxc = locc + int(no_neighbpix_x/2) + 1
  if minr < 0:
    maxr -= minr
    minr = 0
  if maxr > image_xsize:
    minr = minr + (image_xsize - maxr)
    maxr = image_xsize
  if minc < 0:
    maxc -= minc
    minc = 0
  if maxc > image_xsize:
    minc = minc + (image_xsize - maxc)
    maxc = image_xsize  
  neighborpixels_pos = [r*image_xsize+c for r in np.arange(minr, maxr) for c in np.arange(minc, maxc)]
  return neighborpixels_pos, minr, minc

def plotsta(ich, sta_ich, neighborpixels_pos, spkcnt_ich, dt_sta_graphics, time_bins,figpath,spksrc,\
              figformat,maxsta_idx,dt_sta,dur_sta, no_neighbpix_x,locc, locc_, minc, locr, locr_, minr, popt_on,mappedrf_on,\
              minsta_idx, popt_off, mappedrf_off):
    fntsz = 12
    darkyellow = [255./255,200./255,0./255]# [238./255,238./255,0]
    x = np.linspace(0, no_neighbpix_x-1, no_neighbpix_x)
    y = np.linspace(0, no_neighbpix_x-1, no_neighbpix_x)
    x, y = np.meshgrid(x, y) 
    plt.figure(figsize=(8,6))   
    plt.subplot(221)
    plt.plot(sta_ich[:,neighborpixels_pos], color='grey')
    STASHAPES = sta_ich[:,neighborpixels_pos] #sure?
    plt.hold('on')
    #  plt.plot(sta_ich[:,maxsta_idx[1]], color='r')

    #  plt.title('STA - '+ 'spkcnt ' + str(spkcnt_ich) , fontsize=10)
    plt.xlabel('Time (s)', fontsize=fntsz)
    plt.ylabel('Pixel brightness', fontsize=fntsz)
    plt.xticks(np.arange(0, len(sta_ich[:,0]), 1/dt_sta_graphics), time_bins)
    plt.ylim([0,1])
    polish_fig(plt.gca(),xlabcoord=-.15,vis_rightaxis=True, vis_leftaxis=False, vis_topaxis=False, vis_bottomaxis=True, xtick_pos='bottom', ytick_pos='right',xlabel_pos='bottom', ylabel_pos='right',boldvline=False)
    plt.subplots_adjust(hspace=.4)

    # Plotting firing rate*******
    #TODO  
    plt.subplot(222)
    plt.axis('off')

    plt.subplot(223)
    img_on = sta_ich[maxsta_idx[0],neighborpixels_pos]   
    #  sta_f = scipy.io.loadmat('/Users/sahar/Documents/spktanlz/data/retina/P38_06Mar14/ret2/Luiz/neuron_789_maxsta_SA2.mat', squeeze_me=True)
    #  img_on = np.hstack(sta_f['mean_sta'].transpose().reshape([27*27,1]))
    plt.title('Positive peak (t = '+ str(int((maxsta_idx[0]*dt_sta - dur_sta)*1000)) + ' ms)' , fontsize=10, color='k')
    plt.imshow(np.reshape(img_on, [no_neighbpix_x, no_neighbpix_x]), cmap=pylab.cm.gray, interpolation='none')
    POSFIG = np.reshape(img_on, [no_neighbpix_x, no_neighbpix_x])
    POSPEAK = int((maxsta_idx[0]*dt_sta - dur_sta)*1000)
    POSCENTRE = [locc-minc-.5, locr-minr-.5]
    plt.text(locc-minc-.5, locr-minr-.5,'+',color='r',horizontalalignment='center',verticalalignment='center', fontsize=fntsz)

    if popt_on.any():  
        data_fitted_on = Gaussian2D_bivariate((x, y), *popt_on) 
        if data_fitted_on.any():
          plt.contour(x, y, data_fitted_on.reshape(no_neighbpix_x, no_neighbpix_x), [Gaussian2D_bivariate((popt_on[0]+popt_on[3], popt_on[1]+popt_on[4]), *popt_on)], colors='k', linewidths=2)
          plt.text(popt_on[0],popt_on[1],'+',color='k',horizontalalignment='center',verticalalignment='center', fontsize=fntsz)
          if mappedrf_on:
            plt.title('Positive peak (t = '+ str(int((maxsta_idx[0]*dt_sta - dur_sta)*1000)) + ' ms)' , fontsize=fntsz, color=darkyellow)
            plt.contour(x, y, data_fitted_on.reshape(no_neighbpix_x, no_neighbpix_x), [Gaussian2D_bivariate((popt_on[0]+popt_on[3], popt_on[1]+popt_on[4]), *popt_on)], colors='orange', linewidths=2)
            plt.text(popt_on[0],popt_on[1],'+',color='orange',horizontalalignment='center',verticalalignment='center', fontsize=fntsz)
    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.tick_params(axis='y', which='both', right='off', left='off', labelleft='off')

    plt.subplot(224)
    img_off = sta_ich[minsta_idx[0],neighborpixels_pos]   
    plt.title('Negative peak (t = '+ str(int((minsta_idx[0]*dt_sta - dur_sta)*1000)) + ' ms)' , fontsize=fntsz)
    plt.imshow(np.reshape(img_off, [no_neighbpix_x, no_neighbpix_x]), cmap=pylab.cm.gray, interpolation='none')
    NEGFIG = np.reshape(img_off, [no_neighbpix_x, no_neighbpix_x])
    NEGPEAK = int((minsta_idx[0]*dt_sta - dur_sta)*1000)
    NEGCENTRE = [locc-minc-.5, locr-minr-.5]
    plt.text(locc-minc-.5, locr-minr-.5,'+',color='r',horizontalalignment='center',verticalalignment='center', fontsize=fntsz)
    if popt_off.any():
        data_fitted_off = Gaussian2D_bivariate((x, y), *popt_off)
        if data_fitted_off.any():
          plt.contour(x, y, data_fitted_off.reshape(no_neighbpix_x, no_neighbpix_x), [Gaussian2D_bivariate((popt_off[0]+popt_off[3], popt_off[1]+popt_off[4]), *popt_off)], colors='k', linewidths=2)
          plt.text(popt_off[0], popt_off[1],'+',color='k',horizontalalignment='center',verticalalignment='center', fontsize=fntsz)
        #      pdb.set_trace()
          if mappedrf_off:
            plt.title('Negative peak (t = '+ str(int((minsta_idx[0]*dt_sta - dur_sta)*1000)) + ' ms)' , color=tableau[20],fontsize=fntsz)
            plt.contour(x, y, data_fitted_off.reshape(no_neighbpix_x, no_neighbpix_x), [Gaussian2D_bivariate((popt_off[0]+popt_off[3], popt_off[1]+popt_off[4]), *popt_off)], colors='b', linewidths=2)
            plt.text(popt_off[0], popt_off[1],'+',color=tableau[20],horizontalalignment='center',verticalalignment='center', fontsize=fntsz)
    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.tick_params(axis='y', which='both', right='off', left='off', labelleft='off')

    plt.suptitle(spksrc+ str(ich)+ ' - loc (' + "%0.1f" % locr_ + ',' + "%0.1f" % locc_ + ')', fontsize=fntsz, y=.99)    
    figname = os.path.join(figpath,'sta_'+spksrc+str(ich)+'_spkcnt'+str(spkcnt_ich)+'dur'+str(dur_sta)+'dt'+str(dt_sta)+figformat)
    plt.savefig(figname, bbox_inches='tight')  
    print 'figure saved to:', figname    
    plt.close()
    filename = os.path.join(figpath,'sta_'+spksrc+str(ich)+'_spkcnt'+str(spkcnt_ich)+'dur'+str(dur_sta)+'dt'+str(dt_sta)+'.hdf5')
    FL = h5py.File(filename)
    FL.create_dataset('negfig',data=NEGFIG)
    FL.create_dataset('posfig',data=POSFIG)
    FL.create_dataset('negpeak',data=NEGPEAK)
    FL.create_dataset('pospeak',data=POSPEAK)
    FL.create_dataset('poscentre',data=POSCENTRE)
    FL.create_dataset('negcentre',data=NEGCENTRE)
    FL.create_dataset('sta',data=STASHAPES)
    FL.close()

def polish_fig(ax,xlabcoord=-.12,vis_rightaxis=False, vis_leftaxis=True, vis_topaxis=False, vis_bottomaxis=True, xtick_pos='bottom', ytick_pos='left', xlabel_pos='bottom', ylabel_pos='left', boldvline=True, invis_firstlabel=False, a=2):
  mplt.rcParams['xtick.direction'] = 'out'
  mplt.rcParams['ytick.direction'] = 'out'
  mplt.rc('axes',facecolor='ffffff')
  mplt.rc('axes',edgecolor='000000')
  mplt.rc('axes',labelcolor='000000')
  mplt.rc('xtick',color='000000')
  mplt.rc('ytick',color='000000')
  ax.spines["right"].set_visible(vis_rightaxis)
  ax.spines["top"].set_visible(vis_topaxis)
  ax.spines["left"].set_visible(vis_leftaxis)
  ax.spines["bottom"].set_visible(vis_bottomaxis)
  ax.xaxis.set_ticks_position(xtick_pos)
  ax.yaxis.set_ticks_position(ytick_pos)
  ax.xaxis.set_label_position(xlabel_pos)
  ax.yaxis.set_label_position(ylabel_pos)
  ax.xaxis.set_label_coords(0.5, xlabcoord)
      
  if invis_firstlabel:
    plt.setp(ax.get_yticklabels()[0], visible=False)    
    plt.setp(ax.get_xticklabels()[0], visible=False)    

#  ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
#  ax.tick_params(axis='x', pad=15) # distant ticks
  if boldvline:
    ax.axvline(linewidth=a, color="black")
    ax.axhline(linewidth=a, color="black")
    ax.tick_params('both', width=a, which='major', direction='out')


def Gaussian2D_bivariate((x, y), xo, yo, amplitude, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()

# assigning data variables
ch = data['ch']
no_ch = len(np.unique(ch))  
freq = data['freq']
spikes = data['spikes']
loc = data['loc']
loc2d = data['loc2d']
nbins = data['nbins']

# assigning stimulus variables
stimframes = stim['stimframes']
stimpixels = stim['stimpixels']
image_xsize = int(np.sqrt(len(stimpixels[0])))
dur_stim_f = int(stimframes[1]-stimframes[0])
dur_stim = dur_stim_f / float(freq)
print 'The duration of White-noise:', dur_stim*1000, 'ms'

# preparing the binning of the sta array
dur_sta_f = int(dur_sta * freq)
spklatency_sta_f = int(spklatency_sta * freq)
no_bins_sta = int(dur_sta / dt_sta)
no_bins_sta_stim = int(dur_sta / dur_stim)
if no_bins_sta_stim > no_bins_sta:
    raise ValueError('Your required bins for sta (dt_sta) is too small for this stimulus')
no_bins_add = no_bins_sta / no_bins_sta_stim
stimidx_rep = np.repeat(range(len(stimframes)), no_bins_add)
time_bins = np.arange(-dur_sta*10, 1, dt_sta_graphics*10)/10.

stimframes_idxpause = [idx for idx, elm in enumerate(stimframes[:-1]) if (stimframes[idx+1] - elm)>dur_stim_f*2]
print 'frames at which breaks happens (starting from 0):', stimframes_idxpause
spkidx = np.array([], dtype=int)
startf_idx = 0
no_neighbpix_x = np.sqrt(no_neighbpix)
for i, idx_p in enumerate(stimframes_idxpause):
    spkidx_spcf = get_spkidx(data, minamp, stimframes[startf_idx]+dur_sta_f+spklatency_sta_f, stimframes[idx_p]+spklatency_sta_f)
    spkidx = np.append(spkidx, spkidx_spcf)
    startf_idx = idx_p+1
spkidx_spcf = get_spkidx(data, minamp, stimframes[startf_idx]+dur_sta_f+spklatency_sta_f, stimframes[-1]+spklatency_sta_f)
spkidx = np.append(spkidx, spkidx_spcf)

hist_spkcnt = get_hist_spkcnt(data, spkidx)
allspikingChannels = np.where(hist_spkcnt >= minspk)[0] # if all channels spike: ~= np.arange(no_ch)
if slctdchs_sta == 'all':
    spikingChannels = allspikingChannels
else:
    spikingChannels = slctdchs_sta[np.where(hist_spkcnt[slctdchs_sta]>=minspk)[0]]

print 'spiking channels:', spikingChannels

stalatency_on, stalatency_off = (np.zeros(max(spikingChannels)+1), np.zeros(max(spikingChannels)+1))
stapeakresp_on, stapeakresp_off = (np.ones(max(spikingChannels)+1)*.5, np.ones(max(spikingChannels)+1)*.5)
rfdiameters_on, rfdiameters_off = (np.zeros(max(spikingChannels)+1), np.zeros(max(spikingChannels)+1))
starespunits, mappedrfunits = (np.zeros(max(spikingChannels)+1), np.zeros(max(spikingChannels)+1))
popts_on = np.zeros([max(spikingChannels)+1, no_gausfit_params])
popts_off = np.zeros([max(spikingChannels)+1, no_gausfit_params])
locs_spikingChs = np.zeros((max(spikingChannels)+1,2))

# params returned after computation of firing rate  
ffspkcnt2d_avg, ffspkcnt2d_std, ffspkcnt2l_avg, ffspkcnt2l_std = (np.zeros(max(spikingChannels)+1) for i in range(4))
fflatency2dark, fflatency2light, ffduresp2dark, ffduresp2light, ffpeakresp2dark, ffpeakresp2light = (np.zeros(max(spikingChannels)+1) for i in range(6))
ffirstspk2d_avg, ffirstspk2d_std, ffirstspk2l_avg, ffirstspk2l_std = (np.zeros(max(spikingChannels)+1) for i in range(4))
biasidxs = np.ones(max(spikingChannels)+1) * -1.2

for ich in spikingChannels:
    print '\n', spksrc, ':', str(ich), '*************'
    
    # Plot and Compute Bias Index if necessary ******************************
    if calcBiasidx4staFLG:
        pass
        # TODO: its ipython notebook needs to be implemented

    popt_on, popt_off = (np.array([]),np.array([]))
    pcov_on, pcov_off = (np.array([]),np.array([]))
    locr_, locc_ = get_chlocation(ich, spksrc, nbins, loc, coord_spksrc)
    if spksrc == 'cl':
      locr_, locc_ = get_chlocation(ich, spksrc, nbins, loc2d, coord_spksrc) 
    else:
      locr_, locc_ = get_chlocation(ich, spksrc, nbins, loc, coord_spksrc) 
    locs_spikingChs[ich] = np.array([locr_, locc_])

    locr = (locr_ + .5) * image_xsize/64.0 
    locc = (locc_ + .5) * image_xsize/64.0 
    neighborpixels_pos, minr, minc = get_neighborpixels(int(round(locr)), int(round(locc)), image_xsize, no_neighbpix_x)
    spkidx_ich = np.where(ch==ich)[0]
    spkidx_ich = np.intersect1d(spkidx_ich, spkidx)
    sta_ich = np.zeros([no_bins_sta+1, len(stimpixels[0])])

    startfs_sta_ich = spikes[spkidx_ich] - dur_sta_f - spklatency_sta_f
    #joins the paused frames together with the distance dur_stim_f
    startfs_sta_ich_tmp = np.copy(startfs_sta_ich)   
    for pausef_idx in stimframes_idxpause:
      idx_afterp = np.where(startfs_sta_ich > stimframes[pausef_idx])
      startfs_sta_ich_tmp[idx_afterp] -= stimframes[pausef_idx+1] - stimframes[pausef_idx] - dur_stim_f

    startfs_sta_ich_idx = [int(x) for x in (startfs_sta_ich_tmp-stimframes[0])/float(dur_stim_f)*no_bins_add]
    #    pdb.set_trace()    
    for idx, startfi in enumerate(startfs_sta_ich_idx): 
    #      print '\n***\nspike frame:', spikes[spkidx_ich[idx]], '\nstartframe:', startfs_sta_ich_tmp[idx], '\n sta frames:', [stimframes[0]+(i*dur_stim_f/no_bins_add) for i in range(startfi,startfi+no_bins_sta+1)], '\n stimpixels:', stimpixels[stimidx_rep[startfi:startfi+no_bins_sta+1],0],'\n'
      sta_ich += stimpixels[stimidx_rep[startfi:startfi+no_bins_sta+1]]
    sta_ich /= len(spkidx_ich)
    #*************************************************
    # finding positive and negative peak STAs 
    maxsta_idx = np.unravel_index(sta_ich.argmax(), sta_ich.shape)
    minsta_idx = np.unravel_index(sta_ich.argmin(), sta_ich.shape)

    # Calculating the zscores of positive and negative peak STAs: (x-mu)/var  
    zscores_atmaxsta = stats.zscore(sta_ich[maxsta_idx[0],neighborpixels_pos])
    zscores_atminsta = stats.zscore(sta_ich[minsta_idx[0],neighborpixels_pos])
    
    # calculating p-values of peak STAs by computing the cumulative density function at peaks
    # scipy.special.ndtr returns cumulative density function, i.e. the area under the probability density function of the given zscores
    pval_atmaxsta = scipy.special.ndtr(-np.abs(zscores_atmaxsta[maxsta_idx[1]]))
    pval_atminsta = scipy.special.ndtr(-np.abs(zscores_atminsta[minsta_idx[1]]))

    if (checkSTAPeakZFLG and pval_atmaxsta > max_pval_acceptable and pval_atminsta > max_pval_acceptable):
      continue
    #    pdb.set_trace()

    # Evaluate ON response ***************************************************
    mappedrf_on = 0
    x = np.linspace(0, no_neighbpix_x-1, no_neighbpix_x)
    y = np.linspace(0, no_neighbpix_x-1, no_neighbpix_x)
    x, y = np.meshgrid(x, y) 
    img_on = sta_ich[maxsta_idx[0],neighborpixels_pos]   
    # reading sta directly from a file****************
    #    pdb.set_trace()
    #    sta_f = scipy.io.loadmat('/Users/sahar/Documents/spktanlz/data/retina/P38_06Mar14/ret2/Luiz/neuron_789_maxsta_SA2.mat', squeeze_me=True)
    #    img_on = np.hstack(sta_f['mean_sta'].transpose().reshape([27*27,1]))

    initial_guess = [locc-minc-.5, locr-minr-.5, 1, 1, 1, 0, 0]
    #    initial_guess = [locc-minc-1, locr-minr-1, 1, 1, 1, 0, 0]
    try:
      popt_on, pcov_on = opt.curve_fit(Gaussian2D_bivariate, (x,y), img_on, p0=initial_guess)
      popts_on[ich] = popt_on
      popt_onstr = [float("%0.2f" % v) for v in popt_on]
      popt_var_on = [float("%0.2f" % v) for v in pcov_on.diagonal()]
      if PrintGaussParamFLG:
        print '\nON: Gaussian2D was fit by the parameters:\nxo=', popt_onstr[0], '\nyo=', popt_onstr[1], \
              '\namplitude=', popt_onstr[2] , '\nsigma_x=', popt_onstr[3], '\nsigma_y=', popt_onstr[4],'\ntheta=', popt_onstr[5], '\noffset=', popt_onstr[6]     
        print '\nON: Variance of the parameters:\nxo_var=', popt_var_on[0], '\nyo_var=', popt_var_on[1], \
              '\namplitude_var=', popt_var_on[2] , '\nsigma_x_var=', popt_var_on[3], '\nsigma_y_var=', popt_var_on[4],'\ntheta_var=', popt_var_on[5], '\noffset_var=', popt_var_on[6]     

    #      data_fitted_on = Gaussian2D_bivariate((x, y), *popt_on)      
      #pval_atmaxsta < max_pval_acceptable and             
      if (np.abs(popt_var_on[0]) < max_var_acceptable and np.abs(popt_var_on[1]) < max_var_acceptable and popt_onstr[2] > 0):
        print '***On-response RF was mapped ***'
        mappedrf_on = 1
    except:
      print '\n!!!on-response could not be mapped\n'

    # Evaluate OFF response ***************************************************
    mappedrf_off = 0
    img_off = sta_ich[minsta_idx[0],neighborpixels_pos]   
    initial_guess = [locc-minc-.5, locr-minr-.5, 1, 1, 1, 0, 0]

    try:
      popt_off, pcov_off = opt.curve_fit(Gaussian2D_bivariate, (x, y), img_off, p0=initial_guess)
      popts_off[ich] = popt_off
      popt_offstr = [float("%0.2f" % v) for v in popt_off]
      popt_var_off = [float("%0.2f" % v) for v in pcov_off.diagonal()]
      if PrintGaussParamFLG:
        print '\n***\nOFF: Gaussian2D was fit by the parameters:\nxo=', popt_offstr[0], '\nyo=', popt_offstr[1], \
              '\namplitude=', popt_offstr[2] , '\nsigma_x=', popt_offstr[3], '\nsigma_y=', popt_offstr[4],'\ntheta=', popt_offstr[5], '\noffset=', popt_offstr[6] #Gaussian2D_bivariate
        print '\nOFF: Variance of the parameters:\nxo_var=', popt_var_off[0], '\nyo_var=', popt_var_off[1], \
              '\namplitude_var=', popt_var_off[2] , '\nsigma_x_var=', popt_var_off[3], '\nsigma_y_var=', popt_var_off[4],'\ntheta_var=', popt_var_off[5], '\noffset_var=', popt_var_off[6]     
    #      data_fitted_off = Gaussian2D_bivariate((x, y), *popt_off)
      # pval_atminsta < max_pval_acceptable
    #      pdb.set_trace()
      if (np.abs(popt_var_off[0]) < max_var_acceptable and np.abs(popt_var_off[1]) < max_var_acceptable and popt_offstr[2] < 0):
        print '***OFF-response RF was mapped ***'
        mappedrf_off = 1
    except:
      print '\n!!!off-response could not be mapped\n'

    # Setting STA parameters **************************************
    isvalid_tempstaon = (pval_atmaxsta < max_pval_acceptable)
    isvalid_tempstaoff = (pval_atminsta < max_pval_acceptable)
    if isvalid_tempstaon:
      stalatency_on[ich] = (dur_sta - maxsta_idx[0]*dt_sta)*1000
      stapeakresp_on[ich] = sta_ich.max()
    if isvalid_tempstaoff:
      stalatency_off[ich] = (dur_sta - minsta_idx[0]*dt_sta)*1000
      stapeakresp_off[ich] = sta_ich.min()

    #    pdb.set_trace()
    if (isvalid_tempstaon and isvalid_tempstaoff):    
      # if ON peak is nearer to spike than OFF peak
      if maxsta_idx[0] > minsta_idx[0]:
        starespunits[ich] = 1
      else:
        starespunits[ich] = -1
    elif isvalid_tempstaon:
      starespunits[ich] = 1      
    elif isvalid_tempstaoff:
      starespunits[ich] = -1

    # Setting the RF parameters **************************************
    #    pdb.set_trace()
    if mappedrf_on:
      rfdiameters_on[ich] = 2 * np.sqrt(np.abs(popt_on[3] * popt_on[4]))
    if mappedrf_off:
      rfdiameters_off[ich] = 2 * np.sqrt(np.abs(popt_off[3] * popt_off[4]))

    if mappedrf_on and mappedrf_off:    
      # if ON peak is nearer to spike than OFF peak
      if (starespunits[ich] == 1):
        mappedrfunits[ich] = 1
      else:
        mappedrfunits[ich] = -1
    elif mappedrf_on:
      mappedrfunits[ich] = 1
    elif mappedrf_off:
      mappedrfunits[ich] = -1
    
    # receptive field diameters
    rfdiameters_on[ich] *= chip_size * 1000 * 1./image_xsize     
    rfdiameters_off[ich] *= chip_size * 1000 * 1./image_xsize     

    stapeakresp_on[ich] = np.abs(0.5 - stapeakresp_on[ich])
    stapeakresp_off[ich] = np.abs(0.5 - stapeakresp_off[ich])
    #    pdb.set_trace()     
    # Plot STA ***************************************************************
    if staplotFLG and (mappedrf_on or mappedrf_off or not plotJustRespUnitsFLG): 
      plotsta(ich, sta_ich, neighborpixels_pos, hist_spkcnt[ich], dt_sta_graphics, time_bins,figpath, spksrc,\
              figformat,maxsta_idx,dt_sta,dur_sta, no_neighbpix_x,locc, locc_, minc, locr, locr_, minr, popt_on,mappedrf_on,\
              minsta_idx, popt_off, mappedrf_off)
        
    if staSaveParamsFLG:
        print 'saving sta params in a mat file...'
        paramnames = np.array(['Spikingchannels', 'locs_r', 'locs_c', 'biasidxs', 'latency2l', 'latency2d',\
                  'duresp2l', 'duresp2d', 'peakresp2l', 'peakresp2d', 'firstspk2l_avg', 'firstspk2d_avg', \
                  'firstspk2l_std', 'firstspk2d_std', 'spkcnt2l_avg', 'spkcnt2d_avg', 'spkcnt2l_std', 'spkcnt2d_std', \
                  'starespunits', 'stalatency_on', 'stalatency_off', 'stapeakresp_on', 'stapeakresp_off', 'mappedrfunits', 'rfdiameters_on', 'rfdiameters_off', \
                  'rf_on_x0', 'rf_on_y0', 'rf_on_amp', 'rf_on_sigmax', 'rf_on_sigmay', 'rf_on_theta', 'rf_on_offset',\
                  'rf_off_x0', 'rf_off_y0', 'rf_off_amp', 'rf_off_sigmax', 'rf_off_sigmay', 'rf_off_theta', 'rf_off_offset'])
        params = np.transpose(np.vstack([spikingChannels, locs_spikingChs[spikingChannels, 0], locs_spikingChs[spikingChannels, 1], biasidxs[spikingChannels], fflatency2light[spikingChannels], fflatency2dark[spikingChannels],\
                  ffduresp2light[spikingChannels], ffduresp2dark[spikingChannels], ffpeakresp2light[spikingChannels], ffpeakresp2dark[spikingChannels], ffirstspk2l_avg[spikingChannels], ffirstspk2d_avg[spikingChannels], \
                  ffirstspk2l_std[spikingChannels], ffirstspk2d_std[spikingChannels], ffspkcnt2l_avg[spikingChannels], ffspkcnt2d_avg[spikingChannels], ffspkcnt2l_std[spikingChannels], ffspkcnt2d_std[spikingChannels],\
                  starespunits[spikingChannels], stalatency_on[spikingChannels], stalatency_off[spikingChannels], stapeakresp_on[spikingChannels], stapeakresp_off[spikingChannels], mappedrfunits[spikingChannels], rfdiameters_on[spikingChannels], rfdiameters_off[spikingChannels], \
                  popts_on[spikingChannels,0], popts_on[spikingChannels,1],popts_on[spikingChannels,2],popts_on[spikingChannels,3],popts_on[spikingChannels,4],popts_on[spikingChannels,5],popts_on[spikingChannels,6],\
                  popts_off[spikingChannels,0], popts_off[spikingChannels,1],popts_off[spikingChannels,2],popts_off[spikingChannels,3],popts_off[spikingChannels,4],popts_off[spikingChannels,5],popts_off[spikingChannels,6]]))

        params = params.astype(object)
        params = np.insert(params, 0, paramnames, axis=0)
    #    pdb.set_trace()
        paramsf = staparams_matfile+'_'+str(len(spikingChannels))+spksrc+'.mat'
        scipy.io.savemat(paramsf,mdict={'STAfrParams':params}, oned_as='column')
        print 'Parameters were saved into: ', paramsf    


