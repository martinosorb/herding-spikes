import numpy as np
import h5py
import scipy.stats
import scipy.cluster

# This is the functions used in the program that writes the .txt File of the
# online detection into a .hdf File, and marks events detected in multiple
# channels and runs a correlation analysis.


# figure out dimensions of arrays, read spike files
def readSpikesFile(TxtFile, HdfFile, NoisyChFile, recCh=4096,
                   removeCh=0, tMax=0, Sampling=7022):
    if recCh == 4096:
        NCh = 4096
        recCh = np.arange(recCh, dtype=int)
    else:
        NCh = len(recCh)
    b = file(TxtFile + '_online_Spikes.txt')
    X = []
    Amp = []
    Y = []
    for i in b:
        z = np.array(i.split(), dtype=int)
        X.append(z[0])
        Y.append(z[1])
        Amp.append(z[2]/64.)
    b.close()
    SpkCh = np.array(X)
    # want to randomize events that are in same frame.
    SpkT = np.array(Y)+scipy.rand(SpkCh.shape[0])
    SpkAmp = np.array(Amp)
    if tMax == 0:
        tMax = np.ceil(np.max(SpkT)*1./Sampling)
    Ncount = np.histogram(SpkCh, bins=np.arange(NCh+1))[0]
    if removeCh > 0:
        IgnCh = np.nonzero(Ncount[:NCh] > removeCh*tMax)[0]
        Ncount[IgnCh] = 0
    elif removeCh == -1:
        # list of channels to ignore (using the indices used in the spike
        # detection eventually need to be converted)
        g = h5py.File(NoisyChFile, 'r')
        IgnCh = np.array(g['NoisyChannels'].value, dtype=int)
        g.close()
        Ncount[IgnCh] = 0
    else:
        IgnCh = np.array([-1])
    NSpk = np.sum(Ncount)
    PreSelectedEvents = True-np.in1d(SpkCh, IgnCh)
    Ind = np.argsort(SpkT[PreSelectedEvents])  # sort according to time
    Loc = np.zeros((np.sum(PreSelectedEvents), 2), dtype=int)
    Loc[:, 0] = SpkCh[PreSelectedEvents][Ind]/64
    Loc[:, 1] = SpkCh[PreSelectedEvents][Ind] % 64
    g = h5py.File(HdfFile, 'w')
    g.create_dataset('Sampling', data=Sampling)
    g.create_dataset('tMax', data=tMax)
    g.create_dataset('nFrames', data=Sampling*tMax)
    g.create_dataset('recordedChannels', data=recCh)
    if 'RawEvents' in g:
        del g['RawEvents']
    f = g.create_group('RawEvents')
    # to find spikes and raw data in the .txt files
    f.create_dataset('SortInd', data=Ind, dtype=int)
    f.create_dataset('Amplitudes', data=SpkAmp[PreSelectedEvents][Ind])
    f.create_dataset('Channels', data=SpkCh[PreSelectedEvents][Ind])
    f.create_dataset('Times', data=SpkT[PreSelectedEvents][Ind])
    f.create_dataset('Locations', data=Loc)
    f.create_dataset('RepolarizingSpikes',
                     data=np.ones(np.sum(PreSelectedEvents), dtype=bool),
                     dtype=bool)
    g.create_dataset('NoisyChannels', data=IgnCh, dtype=int)
    f.create_dataset('PreSelectedEvents', data=PreSelectedEvents, dtype=bool)
    g.close()
    print(NSpk)
    return


def IsolatedSpikes(HdfFile, DFrames=2, MaxDist=1.5):
    f = h5py.File(HdfFile, 'r+')
    g = f['RawEvents']
    # Sampling = f['Sampling'].value
    Loc = g['Locations'].value
    Amplitudes = g['Amplitudes'].value
    Times = g['Times'].value
    Channels = g['Channels'].value
    RepolarizingSpikes = np.array(g['RepolarizingSpikes'].value, dtype=bool)
    X = np.zeros((1001, 2))
    Xa = np.zeros((1001))
    Xt = np.zeros((1001))
    Ispikes = np.ones(len(Times), dtype=bool)
    for i in range(500):
        X[i, :] = Loc[i, :]
        Xa[i] = Amplitudes[i]
        Xt[i] = Times[i]
    for i in range(500, len(Times)):
        X[i % 1001, :] = Loc[i, :]
        Xa[i % 1001] = Amplitudes[i]
        Xt[i % 1001] = Times[i]
        j = (i-500) % 500
        Ind = np.nonzero(((Xa-Xa[j]) > 0)*(np.abs(Xt-Xt[j]) <= DFrames))[0]
        if any(np.sum((X[Ind, :]-X[j, :][None, :])**2, axis=1) <= MaxDist):
            Ispikes[i] = False
    for i in range(len(Times)-500, len(Times)):
        j = i % 500
        Ind = np.nonzero(((Xa-Xa[j]) > 0)*(np.abs(Xt-Xt[j]) <= DFrames))[0]
        if any(np.sum((X[Ind, :]-X[j, :][None, :])**2, axis=1) <= MaxDist):
            Ispikes[i] = False
    if 'IsolatedSpikes' in g:
        del g['IsolatedSpikes']
    g.create_dataset('IsolatedSpikes', data=Ispikes)
    if 'Locations' in f:
        del f['Locations']
    f.create_dataset('Locations', data=Loc[Ispikes, :]+0.5)
    if 'Amplitudes' in f:
        del f['Amplitudes']
    f.create_dataset('Amplitudes', data=Amplitudes[Ispikes])
    if 'Times' in f:
        del f['Times']
    f.create_dataset('Times', data=Times[Ispikes])
    if 'Channels' in f:
        del f['Channels']
    f.create_dataset('Channels', data=Channels[Ispikes])
    if 'RepolarizingSpikes' in f:
        del f['RepolarizingSpikes']
    f.create_dataset('RepolarizingSpikes', data=RepolarizingSpikes[Ispikes])
    if 'IncludeLongSpikes' in f:
        del f['IncludeLongSpikes']
    f.create_dataset('IncludeLongSpikes', data=False, dtype=bool)
    f.close()
    return


def CorrelationAnalysis1(HdfFile, NextN=100, NPoissonNoise=120,
                         dClockspikes=4, Nignore=20):
    f = h5py.File(HdfFile, 'r')
    g = f['RawEvents']
    Sampling = f['Sampling'].value
    Loc = f['Locations'].value
    tMax = f['tMax'].value
    # nFrames = f['nFrames'].value
    Times = f['Times'].value
    Spikes = f['Channels'].value
    SingleSpk = np.array(f['RepolarizingSpikes'].value, dtype=bool)
    f.close()

    H0 = np.histogram(Spikes, bins=np.arange(64**2+1))[0]
    CountNew = H0.copy()
    # CountNew[1:]+=H0[:-1]
    # CountNew[128:]+=H0[:-128]
    # CountNew[129:]+=H0[:-129]
    # find local maxima---------------------------------
    ImgRes = 64
    SHistAvg = np.histogram2d(Loc[SingleSpk, 0], Loc[SingleSpk, 1],
                              bins=(np.arange(65), np.arange(65)))[0]
    SHistAvgMax = np.zeros(SHistAvg.shape, dtype=int)
    for i in range(2):
        for j in range(i, 2):
            if np.sqrt(i**2+j**2) < 1.5:
                if j == 0:
                    print(j)
                elif i == j:
                    SHistAvgMax[i:, i:] += SHistAvg[:-i, :-i] < \
                        SHistAvg[i:, i:]
                    SHistAvgMax[:-i, :-i] += SHistAvg[i:, i:] < \
                        SHistAvg[:-i, :-i]
                    SHistAvgMax[i:, :-i] += SHistAvg[:-i, i:] < \
                        SHistAvg[i:, :-i]
                    SHistAvgMax[:-i, i:] += SHistAvg[i:, :-i] < \
                        SHistAvg[:-i, i:]
                elif i == 0:
                    SHistAvgMax[:, j:] += SHistAvg[:, :-j] < SHistAvg[:, j:]
                    SHistAvgMax[j:, :] += SHistAvg[:-j, :] < SHistAvg[j:, :]
                    SHistAvgMax[:, :-j] += SHistAvg[:, j:] < SHistAvg[:, :-j]
                    SHistAvgMax[:-j, :] += SHistAvg[j:, :] < SHistAvg[:-j, :]
    SHistAvgMax = (SHistAvgMax).flatten()
    XX = np.array(np.nonzero((SHistAvgMax >= 3)*(SHistAvg.flatten() >
                  np.percentile(SHistAvg, 10)))[0], dtype=int)
    Freq = SHistAvg.flatten()[XX]
    Ncomp = len(XX)-Nignore
    print(Ncomp)
    XX = XX[np.argsort(np.argsort(
            np.abs(np.log(Freq*1./np.mean(Freq))))) < Ncomp]
    SHistAvgMax = np.zeros(64**2, dtype=int)
    for i in np.array([0]):
        SHistAvgMax[XX+i] = np.arange(1, len(XX+1))
    # assign local maxima to units
    Units = SHistAvgMax[np.array(Loc[:, 1], dtype=int) +
                        np.array(Loc[:, 0], dtype=int)*64]
    # count spikes
    # NSpikes0R=np.histogram(Units,bins=np.arange(Ncomp+2))[0]
    NChannels = len(CountNew)  # +1
    AddNoise = int(np.clip((tMax*NPoissonNoise)/3600, 1, 1e12))
    NSpikes = CountNew + AddNoise
    MeanS = np.mean(NSpikes)
    # NTotal=int(np.sum(CountNew))#int(np.sum(NSpikes[1:]-AddNoise))
    ClockSpikes = np.arange(0, tMax*Sampling, dClockspikes)
    Ncs = int((tMax*Sampling-1)/dClockspikes)+1
    # need those for later (telling where Poisson and real spikes should go for
    # next analysis, no need for amplitudes)
    # VSpikes=np.concatenate((np.ones(Spikes.shape[0],dtype=bool),
    #         np.zeros(NChannels*AddNoise,dtype=bool)))
    PSpikes = np.concatenate((np.zeros(Spikes.shape[0], dtype=bool),
                             np.ones(NChannels*AddNoise, dtype=bool)))
    Spikes = np.concatenate(
            (Spikes,
             np.repeat(np.arange(NChannels), AddNoise)[np.argsort(
                scipy.rand(AddNoise * (NChannels)))],
                np.zeros(Ncs, dtype=int)))
    Units = np.concatenate((Units, np.zeros(NChannels*AddNoise, dtype=int),
                            np.zeros(Ncs, dtype=int)-1))
    # SpikesT=np.concatenate((Times,\
    # np.array(scipy.rand(AddNoise*(NChannels))*Sampling*tMax,dtype=int),
    # ClockSpikes))
    print(len(Times), Ncs, len(ClockSpikes), AddNoise)
    SpikesT = np.concatenate(
        (Times, np.interp(scipy.rand(AddNoise*(NChannels))*(len(Times)+Ncs-1),
         np.arange(len(Times)+Ncs),
         np.sort(np.concatenate((Times, ClockSpikes)))), ClockSpikes))
    SortedT = np.argsort(SpikesT)
    Spikes = Spikes[SortedT]
    Units = Units[SortedT]
    SpikesT = np.sort(SpikesT)
    # TotSpikes=np.sum(NSpikes0R)
    # events with recalibrations for correlation index
    UnitsXI = np.nonzero(Units)[0]  # exclude Poissonspikes
    UnitsX = np.clip(Units[UnitsXI], 0, len(XX))
    # Correlation matrix
    X = np.zeros((NChannels, Ncomp+1), dtype=int)
    # where there is no reference spike, how many reference spikes (histogram)
    NZero = np.zeros((NChannels))
    k = NextN-1
    for i in range(UnitsXI[NextN], UnitsXI[-NextN]):
        if UnitsXI[k] < i:
            k += 1
            a = np.unique(UnitsX[k-NextN:k+NextN+1])
            b = 2*NextN-len(a)+1  # np.sum(Units[i-NextN:i+NextN+1]==0)
        NZero[Spikes[i]] += b
        X[Spikes[i], a] += 1
    NSpikesR = np.sum(X, axis=0)
    NSpikesR[0] = 0
    TotSpikes = np.sum(NSpikesR[1:])
    MeanSR = np.mean(NSpikesR[1:])
    SpikeHisto = np.zeros((NChannels, Ncomp+1))
    for i in range(1, NChannels):
        SpikeHisto[i, :] = NSpikes[i]*NSpikesR*NextN*2. / \
            (TotSpikes-(NSpikes[i]*(SHistAvgMax[i] != 0)))  # ok
    # remove connections between units with different activities
    Z = np.zeros((NChannels, Ncomp+1), dtype=bool)
    for i in range(1, NChannels):
        SpikeHistoTest = (SpikeHisto[i, :]*(1.-(NZero[i]/2./NextN /
                          (NSpikes[i]+(NSpikes[i] == 0)))))
        SpikeHistoTest[SpikeHistoTest <= 0.01] = 0.01
        Z[i, :] = (scipy.stats.poisson.sf(X[i, :]+0.5, SpikeHistoTest) <
                   0.999999)
    Z[0, :] = False
    Z[:, 0] = False
    print(np.sum(Z)/(NChannels-1.))
    SpikeHisto = np.zeros((NChannels, Ncomp+1))
    for i in range(1, NChannels):
        fracGood = (np.sum(NSpikesR[Z[i, :]])+NSpikesR[0])*1./np.sum(NSpikesR)
        SpikeHisto[i, :] = NSpikes[i]*NSpikesR*NextN*2. / \
            (TotSpikes-(NSpikes[i]*(SHistAvgMax[i] != 0)))/fracGood
    # units with increasing # of connections
    NCorr = np.zeros((NChannels, Ncomp+1), dtype=bool)
    probC = np.array([0.001])
    a = XX.copy()  # np.arange(1,1000,dtype=int)
    for k in range(len(probC)):
        for i in range(1, NChannels):
            b = np.arange(1, Ncomp+1)[(
                ((a/ImgRes-i/ImgRes)**2 + (a % ImgRes-i % ImgRes)**2) <
                ImgRes**2/2)*(((a/ImgRes-i/ImgRes)**2 +
                              (a % ImgRes-i % ImgRes)**2) > 1)]
            SpikeHistoTest = (SpikeHisto[i, b]*(1.-(NZero[i]/2./NextN /
                              (NSpikes[i]+(NSpikes[i] == 0)))))
            SpikeHistoTest[SpikeHistoTest <= 0.001] = 0.001
            A = (scipy.stats.poisson.sf(X[i, b]-0.5, SpikeHistoTest) <
                 (probC[k]*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1. /
                                   np.clip(NSpikesR[b], 1, 1e12), 1e-12, 10))))
            B = (scipy.stats.poisson.sf(X[i, b]+0.5, SpikeHistoTest) <
                 (probC[k]*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1. /
                                   np.clip(NSpikesR[b],
                                   1, 1e12), 1e-12, 10))))*(X[i, b] > 0)
            Q = (scipy.stats.poisson.sf(X[i, b]-0.5, SpikeHistoTest) -
                 scipy.stats.poisson.sf(X[i, b]+0.5, SpikeHistoTest))[B-A]
            A[B-A] = scipy.rand(np.sum(B-A)) < (Q/2.)
            NCorr[i, b] = A
        NCorr[0, :] = False
        NCorr[:, 0] = False
        print(np.sum(NCorr)/(NChannels))
    ###
    NCorrClusterLinkage = scipy.cluster.hierarchy.ward(NCorr)
    NCorrFCluster = scipy.cluster.hierarchy.fcluster(NCorrClusterLinkage,
                                                     11, 'maxclust')
    NCorrX = np.zeros((NChannels, Ncomp+1), dtype=bool)
    for k in range(len(probC)):
        for i in range(1, NChannels):
            b = np.arange(1, Ncomp+1)[(((a/ImgRes-i/ImgRes) ** 2 +
                                      (a % ImgRes-i % ImgRes) ** 2) <
                                      ImgRes**2/16)*(((a/ImgRes-i/ImgRes)**2 +
                                                     (a % ImgRes-i %
                                                      ImgRes)**2) > 1)]
            SpikeHistoTest = (SpikeHisto[i, b]*(1.-(NZero[i]/2./NextN /
                              (NSpikes[i]+(NSpikes[i] == 0)))))
            SpikeHistoTest[SpikeHistoTest <= 0.001] = 0.001
            A = (scipy.stats.poisson.sf(X[i, b]-0.5, SpikeHistoTest) <
                 (0.25*probC[k]*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1. /
                  np.clip(NSpikesR[b], 1, 1e12), 1e-12, 10))))
            B = (scipy.stats.poisson.sf(X[i, b]+0.5, SpikeHistoTest) <
                 (0.25*probC[k]*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1. /
                  np.clip(NSpikesR[b], 1, 1e12), 1e-12, 10))))*(X[i, b] > 0)
            Q = (scipy.stats.poisson.sf(X[i, b]-0.5, SpikeHistoTest) -
                 scipy.stats.poisson.sf(X[i, b]+0.5, SpikeHistoTest))[B-A]
            A[B-A] = scipy.rand(np.sum(B-A)) < (Q/2.)
            NCorrX[i, b] = A
    NCorrX[0, :] = False
    NCorrX[:, 0] = False
    print(np.sum(NCorrX)/(NChannels))
    ###
    C = np.zeros(len(Spikes))
    XY = np.concatenate((np.array([-1]), XX))
    k = NextN-1
    for i in range(UnitsXI[NextN], UnitsXI[-NextN]):
        if UnitsXI[k] < i:
            k += 1
            x0 = np.array(UnitsX[k-NextN:k+NextN+1], dtype=int)
            y = np.array(Spikes[UnitsXI[k-NextN:k+NextN+1]], dtype=int)
        if Spikes[i] != 0:
            a = int(Spikes[i])
            x = x0*Z[a, x0]
            # x=np.array(Units[i-NextN:i+NextN+1],dtype=int)
            # x=x*Z[a,x]
            # y=np.array(Spikes[i-NextN:i+NextN+1],dtype=int)#
            a1 = np.array(np.unique(x), dtype=int)[1:]
            a2 = np.array(np.unique(x[NextN/2:-NextN/2+1]), dtype=int)[1:]
            a3 = a1[(((a/ImgRes-XY[a1]/ImgRes)**2 +
                    (a % ImgRes-XY[a1] % ImgRes)**2) < 64**2/16) *
                    (((a/ImgRes-XY[a1]/ImgRes)**2 +
                      (a % ImgRes-XY[a1] % ImgRes)**2) > 1)]
            a4 = a2[(((a/ImgRes-XY[a2]/ImgRes)**2 +
                    (a % ImgRes-XY[a2] % ImgRes)**2) < 64**2/16) *
                    (((a/ImgRes-XY[a2]/ImgRes)**2 +
                      (a % ImgRes-XY[a2] % ImgRes)**2) > 1)]
            a1 = a1[(((a/ImgRes-XY[a1]/ImgRes)**2 +
                    (a % ImgRes-XY[a1] % ImgRes)**2) < 64**2/2) *
                    (((a/ImgRes-XY[a1]/ImgRes)**2 +
                      (a % ImgRes-XY[a1] % ImgRes)**2) > 1)]
            a2 = a2[(((a/ImgRes-XY[a2]/ImgRes)**2 +
                    (a % ImgRes-XY[a2] % ImgRes)**2) < 64**2/2) *
                    (((a/ImgRes-XY[a2]/ImgRes)**2 +
                      (a % ImgRes-XY[a2] % ImgRes)**2) > 1)]
            a5 = np.unique(y[np.in1d(x, a1)])
            a6 = np.unique(y[np.in1d(x, a2)])
            a7 = np.unique(y[np.in1d(x, a3)])
            a8 = np.unique(y[np.in1d(x, a4)])
            c1 = 0
            c2 = 0
            c3 = 0
            c4 = 0
            if len(a4) > 4:
                b = a4[NCorrX[a, a4]]
                c = a8
                if len(b) >= 5:
                    c4 = (np.sum(NCorrX[c, :][:, b]) -
                          len(c))*(1./len(a4))/np.clip(len(c), 1, 1e6)
            if len(a3) > 4:
                b = a3[NCorrX[a, a3]]
                c = a7
                if len(b) >= 5:
                    c3 = (np.sum(NCorrX[c, :][:, b]) -
                          len(c))*(1./len(a3))/np.clip(len(c), 1, 1e6)
            if len(a2) > 4:
                b = a2[NCorr[a, a2]]
                c = a6
                if len(b) >= 4:
                    c2 = (np.sum(NCorr[c, :][:, b]) -
                          len(c))*(1./len(a2))/np.clip(len(c), 1, 1e6)
            if len(a1) > 3:
                b = a1[NCorr[a, a1]]
                c = a5
                if len(b) >= 4:
                    c1 = (np.sum(NCorr[c, :][:, b]) -
                          len(c))*(1./len(a1))/np.clip(len(c), 1, 1e6)
            C[i] = max(C[i], c1, c2, c3, c4)
    g = h5py.File(HdfFile, 'r+')
    SpikesS = np.zeros(Spikes.shape, dtype=int)
    SpikesCS = np.zeros(Spikes.shape)
    SpikesTS = np.zeros(Spikes.shape)
    SpikesTS[SortedT] = SpikesT
    SpikesCS[SortedT] = C
    SpikesS[SortedT] = Spikes
    L = PSpikes.shape[0]
    if 'CorrelationAnalysis' in g:
        del g['CorrelationAnalysis']
    if 'Units' in g:
        del g['Units']
    if 'Corr' in g:
        del g['Corr']
    f = g.create_group('CorrelationAnalysis')
    f.create_dataset('Cluster', data=NCorrFCluster)
    g.create_dataset('Units', data=SpikesS[:L][True-PSpikes])
    f.create_dataset('AllUnits', data=SpikesS[:L])
    f.create_dataset('PoissonUnits', data=SpikesS[:L][PSpikes])
    f.create_dataset('PoissonTimes', data=SpikesTS[:L][PSpikes])
    g.create_dataset('Corr', data=SpikesCS[:L][True-PSpikes])
    f.create_dataset('PoissonCorr', data=SpikesCS[:L][PSpikes])
    f.create_dataset('PoissonInd', data=PSpikes)
    if 'Parameter' in f:
        del f['Parameter']
    h = f.create_group('Parameter')
    h.create_dataset('LocMax', data=XX)
    h.create_dataset('MaximaResolution', data=ImgRes)
    h.create_dataset('NextN', data=NextN)
    h.create_dataset('NPoissonNoise', data=NPoissonNoise)
    h.create_dataset('dClockspikes', data=dClockspikes)
    h.create_dataset('Nignore', data=Nignore)
    g.close()
    return


def CorrelationAnalysis2(HdfFile):
    f = h5py.File(HdfFile, 'r')
    g = f['CorrelationAnalysis']
    # Sampling = f['Sampling'].value
    # tMax = f['tMax'].value
    # nFrames = f['nFrames'].value
    CS = f['Corr'].value
    CPoisson = g['PoissonCorr'].value
    PoissonInd = g['PoissonInd'].value
    SpikesS = f['Units'].value
    SpikesPoisson = g['PoissonUnits'].value
    SpikesAll = g['AllUnits'].value
    # SpikesTS = f['Times'].value
    SpikesAmpS = np.clip(f['Amplitudes'].value, 0.1, 1e6)
    SingleSpk = np.array(f['RepolarizingSpikes'].value, dtype=bool)
    IncludeLongSpikes = np.array(f['IncludeLongSpikes'].value, dtype=bool)
    f.close()
    LIndS = np.zeros((SpikesS.shape[0], 1+IncludeLongSpikes), dtype=bool)
    LInd = np.zeros((SpikesAll.shape[0], 1+IncludeLongSpikes), dtype=bool)
    LIndS[:, 0] = SingleSpk  # short spikes
    LInd[True-PoissonInd, 0] = SingleSpk  # short spikes
    if IncludeLongSpikes:
        LIndS[:, 1] = (True-SingleSpk)  # long spikes
        LInd[True-PoissonInd, 1] = (True-SingleSpk)  # long spikes
    NChannels = 64**2  # np.sum(MaskX<>0)+1
    Pnew = np.zeros(len(CS.flatten()))
    NSpikes = np.histogram2d(SpikesS,
                             (True-SingleSpk), bins=(np.arange(NChannels+1),
                                                     np.arange(3)))[0]
    AmpCiKink = np.zeros((NChannels, 1+IncludeLongSpikes))
    Pval = np.zeros((NChannels, 1+IncludeLongSpikes))+1.
    CINoise = np.zeros((NChannels, 1+IncludeLongSpikes))+1.
    fNoise = np.zeros((NChannels, 1+IncludeLongSpikes))+1.
    for ii in range(1+IncludeLongSpikes):
        SameSpikes = (PoissonInd-LInd[:, ii]) != 0
        for i in range(NChannels):
            if (NSpikes[i, ii]) > 50:
                PoissonIndx = PoissonInd[SpikesAll == i][
                               SameSpikes[SpikesAll == i]]
                iInd = True-PoissonIndx
                ChAmp = np.zeros(iInd.shape[0])
                ChC = np.zeros(iInd.shape[0])
                ChAmp[iInd] = SpikesAmpS[(SpikesS == i)*LIndS[:, ii]]
                ChC[iInd] = CS[(SpikesS == i)*LIndS[:, ii]]
                ChC[True-iInd] = CPoisson[(SpikesPoisson == i)]
                # sort CI for Poisson spikes
                N = ChC.shape[0]
                cP = CPoisson[SpikesPoisson == i].flatten() + 1e-6*scipy.rand(
                    np.sum(SpikesPoisson == i))
                aP = np.argsort(cP)
                # sort CI for putative events (where spikes are in
                # distribution of Poisson events)
                Rnd = 1e-6*scipy.rand(N)
                aX = np.argsort(ChC.flatten()+Rnd)
                aXS = np.argsort(np.argsort(aX[True-PoissonIndx]))
                iX = np.interp(ChC.flatten()+Rnd,
                               np.concatenate((np.array([0]), cP[aP])),
                               np.linspace(0, 1, len(aP)+1))
                iXS = iX[True-PoissonIndx]  # values between [0,1]
                # sort amplitudes for putative events
                Rnd = 1e-6*scipy.rand(N)
                aY = np.argsort(ChAmp.flatten()+Rnd)
                aYS = np.argsort(np.argsort(aY[True-PoissonIndx]))
                rY = np.argsort(np.argsort(ChAmp.flatten()+Rnd))
                # cumulative sums for KStest (is that one allowed?)
                X = 1.*np.cumsum(rY[aX])
                X *= 0.5/np.clip(X[-1], 1, 1e12)
                Y = np.cumsum(iX[aY])
                Y *= 0.5/np.clip(Y[-1], 1, 1e12)
                NS = np.sum((SpikesS == i)*(LIndS[:, ii]))
                # Ncount[i]=NS
                # AmpCiArea[i]=np.min(X+Y-(np.arange(1,N+1)*1./N))
                # map for excluding Poisson spikes???
                PmapX = np.arange(N)-np.cumsum(
                    0.5*((ChAmp[aX] <= 0.05)+(ChAmp[aY] <= 0.05)))
                BiasC = np.clip(np.sqrt(np.linspace(2./N, 2, N)*np.linspace(
                    2, 2./N, N)), 0.1, 1)
                AmpCiKink[i, ii] = np.mean(PmapX[np.argsort(
                    (X+Y-np.arange(1, N+1)*1./N)*1./BiasC)][:NS/20])*1./(NS-1)
                CINoise[i, ii] = 1.-np.mean(
                    np.clip((np.sort(iXS)*NS-np.arange(NS)) *
                            1./(np.arange(1, NS+1)[::-1]), 0, 1))
                # KStest
                Pval[i, ii] = scipy.stats.kstest(X+Y, 'uniform',
                                                 alternative='greater',
                                                 mode='approx')[1]
                if Pval[i, ii] < 0.01:
                    fNoise[i, ii] = AmpCiKink[i, ii]
                    print(i, Pval[i, ii])
                else:
                    fNoise[i, ii] = CINoise[i, ii]
                # have to decide which spikes are noisy
                # know the fraction of noise and the distribution of noise
                # need to average; use reflecting boundaries
                NPoisson = len(aP)
                # have (1-fNoise) space to increase Poisson probability,
                # or as factor, mean(pPoisson)/fNoise
                # should be bound on pPoisson
                # fPoisson = NPoisson*1./NS # COMMENTED LATER
                nNextC = (NPoisson+NS)
                Cp = np.argsort(np.concatenate(
                    (CPoisson[SpikesPoisson == i].flatten(),
                     CS[(SpikesS == i)*LIndS[:, ii]].flatten()))) < NPoisson
                # changed def
                Cpc = np.cumsum(np.concatenate(
                    (np.array([0]),
                     Cp[:nNextC][::-1], Cp, Cp[-nNextC:][::-1])))
                Ind = np.nonzero(Cp == 0)[0]
                jj = 0.025
                w = True
                while jj <= 1. and w:
                    jj += 0.025
                    nAll = int(nNextC*jj)
                    pPoisson = (Cpc[Ind+nNextC+nAll+1] -
                                Cpc[Ind+nNextC-nAll])/(2.*nAll+1.)+1e-6
                    # should allow for being 10% off
                    w = (np.max(pPoisson)*fNoise[i, ii]) > 1.1*np.mean(
                                                                      pPoisson)
                nS = int(jj*NS)
                pPoisson = (Cpc[Ind+nNextC+nAll+1] -
                            Cpc[Ind+nNextC-nAll])/(2.*nAll+1.)+1e-6
                # print np.mean(pPoisson), fPoisson
                # does not go well with clipping...
                pPoisson = np.clip(
                    pPoisson*fNoise[i, ii]*1./np.mean(pPoisson), 0, 1)
                # do a round of adjustment to fNoise
                pPoisson = np.clip(
                     pPoisson*fNoise[i, ii]*1./np.mean(pPoisson), 0, 1)
                # do a round of adjustment to fNoise
                pPoisson = np.clip(
                     pPoisson*fNoise[i, ii]*1./np.mean(pPoisson), 0, 1)
                fTP0 = np.zeros(NS)
                fTP0[aXS] = 1.-pPoisson
                # now have to do the same with amplitudes
                # Amplitude_CI dependence
                AmpRefl = np.cumsum(np.concatenate(
                    (np.array([0]),
                     fTP0[aYS][:nS][::-1], fTP0[aYS], fTP0[aYS][-nS:][::-1])))
                # normalize to one
                pAmp = (AmpRefl[(2*nS+1):]-AmpRefl[:-(2*nS+1)])/(2.*nS+1.)+1e-6
                # print np.mean(pAmp)
                pAmp *= 1./np.mean(pAmp)  # does not go well with clipping...
                # but shouldn't have that a large effect
                # modify probabilities multiplicatively and locally
                Cpc2 = np.cumsum(np.concatenate(
                    (np.array([0]), fTP0[:nS][::-1], fTP0, fTP0[-nS:][::-1])))
                fTP1 = np.zeros(NS)
                # will correct for clipping in the normalization step
                fTP1[aYS] = np.clip(fTP0[aYS]*pAmp, 0, 1)
                Cpc3 = np.cumsum(np.concatenate(
                    (np.array([0]), fTP1[:nS][::-1], fTP1, fTP1[-nS:][::-1])))
                Cnew = (Cpc3[2*nS+1:]-Cpc3[:-2*nS-1])*0.5/(nS+0.5)
                Cold = (Cpc2[2*nS+1:]-Cpc2[:-2*nS-1])*0.5/(nS+0.5)
                # correct for a local amplitude bias (in CI)
                fTP1 *= Cold*1./np.clip(Cnew, 1e-12, 1e12)
                Pnew[(SpikesS == i)*LIndS[:, ii]] = np.clip(fTP1, 0, 1)
    # save results
    g = h5py.File(HdfFile, 'r+')
    f = g['CorrelationAnalysis']
    if 'Probability' in f:
        del f['Probability']
    if 'Noise' in f:
        del f['Noise']
    if 'pValue' in f:
        del f['pValue']
    f.create_dataset('Probability', data=np.array(Pnew))
    f.create_dataset('Noise', data=np.array(fNoise))
    f.create_dataset('pValue', data=np.array(Pval))
    g.close()
    return
