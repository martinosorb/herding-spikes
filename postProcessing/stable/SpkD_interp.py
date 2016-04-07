import numpy as np
import h5py
import scipy.stats
import scipy.cluster
import time

###This is the functions used in the program that writes the .txt Files into a .hdf Files
# and estimates the spatial origins of spikes, marks events detected in multiple channels,
# runs a correlation analysis and, based on these results, does a clustering.

### read Info file
def readInfoFile(TxtFile, HdfFile):
	X=[]
	Y=[]
	recCh=[]
	k=0
	n=0
	n1=0
	b=file(TxtFile + '_Info.txt')
	for i in b:
		if '#' in i:
			k+=1
			print i
			if k==8:
				XCh4=np.zeros((len(recCh),12),dtype=int)
				n0=0
			if k==9:
				XCh5=np.zeros((len(recCh),9),dtype=int)
				SIprod=np.zeros((len(recCh),13),dtype=long)
				QdAvg=np.zeros((len(recCh)),dtype=int)
				n0=0
		else:
			if k==1:
				nFrames=int(i)
				print i
			if k==2:
				tMax=float(i)
				print i
			if k==3:
				Sampling=float(i[:-1])
				print i
			if k==4:
				SpkThreshold=np.array(i.split(),dtype=int)
				print i
			if k==5:
				SpkRPthreshold=int(i)
				print i
			if k==6:
				recalibrationTrigger=int(i)
				print i
			if k==7:
				recCh.append(int(i))
			if k==8:
				XCh4[n0,:]=np.array(i.split(),dtype=int)
				n0+=1
			if k==9:
				XCh5[n0,:]=np.array(i.split(),dtype=int)
				n0+=1
			if k==10:
				z=np.array(i.split(),dtype=int)
				X.append(z[0])
				QdAvg+=z[1:]
				n+=1
			if k==11:
				SqIglobal=np.array(i.split(),dtype=long)
			if k==12:
				SqIv=np.array(i.split(),dtype=long)
			if k==13:
				SIprod[:,n1]=np.array(i.split(),dtype=long)
				n1+=1
			if k==14:
				Vsbias=-np.array(i.split(),dtype=long)/64./nFrames#in ADCcounts
			if k==15:
				Vsqbias=np.array(i.split(),dtype=long)/64.**2/nFrames
	b.close()
	Recalib=np.array(X)# was used for recalibration events
	Qdavg=QdAvg/64./n#in ADCcounts (for individual channels!)
	recCh=np.array(recCh,dtype=int)
	NCh=len(recCh)
	f=h5py.File(HdfFile,'w')
	#Parameters of the recording and spike detection
	f.create_dataset('Sampling', data=Sampling)
	f.create_dataset('tMax', data=tMax)
	f.create_dataset('nFrames', data=nFrames)
	f.create_dataset('AmplitudeThresholds', data=SpkThreshold)
	f.create_dataset('RepolarizationThreshold', data=SpkRPthreshold)
	f.create_dataset('ignoreRecalibrations', data=recalibrationTrigger)
	#f.create_dataset('lenCutouts', data=Ncut)
	#f.create_dataset('lenCutoutsLong', data=NcutL)
	#f.create_dataset('PeakIndex', data=Npeak)
	#used channels, noise and recalibration events
	#----------------------------------------------------------
	f.create_dataset('recordedChannels', data=recCh)#need this to eventually convert a subset of channels into 64x64 grid
	f.create_dataset('RecalibrationEvents', data=Recalib)# timestamps of recalibration events
	f.create_dataset('ChannelVariability', data=Qdavg)# variability estimate of the channel noise for each channel
	#f.create_dataset('NoisyChannels', data=IgnCh)#channels which I have ignored in the analysis
	#global fluctuations in voltage
	#----------------------------------------------------
	g=f.create_group('GlobalVoltageFluctuations')
	#g.create_dataset('medianVoltage', data=B[:-20]/3.)
	#g.create_dataset('medianVoltageStd', data=Bstd)
	g.create_dataset('sumSqVglobal', data=SqIglobal)#sum of squared increments of median voltage
	g.create_dataset('sumSqVchannels', data=SqIv)#sum of squared increments of individual voltage traces
	g.create_dataset('sumVproduct', data=SIprod)#sum of the product of both
	g.create_dataset('VFBias', data=Vsbias)#How much a linear predictor would take into account global fluctuations
	#i.e. fraction of global fluctuations that would explain (signal-global fluctuations) best
	g.create_dataset('StdVFBias', data=np.sqrt(Vsqbias))#standard deviation of this
	f.close()
	return NCh, tMax, Sampling
	
	
### Bias correction (median voltage)
def readAvgFile(TxtFile, HdfFile):
	b=file(TxtFile + '_Avg.txt')
	X=[]
	for i in b:
		X.append(int(i))
	b.close()
	B=np.array(X)
	Bstd=np.std(B)/3.
	f=h5py.File(HdfFile,'r+')
	g=f['GlobalVoltageFluctuations']
	g.create_dataset('medianVoltage', data=B/3.)
	g.create_dataset('medianVoltageStd', data=Bstd)
	f.close()
	return
	
### figure out dimensions of arrays, read spike files
def readSpikesFile(TxtFile, HdfFile, NoisyChFile, NCh, removeCh, tMax):
	b=file(TxtFile + '_Spikes.txt')
	X=[]
	Amp=[]
	Y=[]
	Y1=[]
	Y2=[]
	YR=[]
	for i in b:
		z=np.array(i.split(),dtype=int)
		X.append(z[0])
		Y.append(z[1])
		Amp.append((z[2])*1./z[3])
		Y1.append(z[4])
		YR.append(z[5])
		Y2.append(9+(np.sum((z[6:9]+z[7:])/2)+(z[6]+z[9])/2)*3)
	b.close()
	b=file(TxtFile + '_SpikesX.txt')
	for i in b:
		z=np.array(i.split(),dtype=int)
		X.append(z[0]+NCh)
		Y.append(z[1])
		Amp.append((z[2])*1./z[3])
		Y1.append(z[4])
		YR.append(z[5])
		Y2.append(12+(np.sum(z[5:9])))
	b.close()
	SpkCh=np.array(X)
	SpkT=np.array(Y)+scipy.rand(SpkCh.shape[0])#want to randomize events that are in same frame.
	SpkAmp=np.array(Amp)
	PeakL=np.array(Y1,dtype=bool)
	Ntraces=np.array(Y2)#not using that one here, maybe should store it anyway
	RecalibOffset=np.array(YR)
	Ncount=np.histogram(SpkCh,bins=np.arange(2*NCh+1))[0]
	if removeCh>0:
		IgnCh=np.nonzero(NCount[:NCh]>removeCh*tMax)[0]
		Ncount[IgnCh]=0
	elif removeCh==-1:
		###list of channels to ignore (using the indices used in the spike detection eventually need to be converted)
		g=h5py.File(NoisyChFile,'r')
		IgnCh=np.array(g['NoisyChannels'].value,dtype=int)
		g.close()
		Ncount[IgnCh]=0
	else :
		IgnCh=np.array([-1])
	NSpk=np.sum(Ncount)
	PreSelectedEvents=True-np.in1d(SpkCh,IgnCh)
	Ind=np.argsort(SpkT[PreSelectedEvents])#sort according to time
	g=h5py.File(HdfFile,'r+')
	if 'RawEvents' in g:
		del g['RawEvents']
	f=g.create_group('RawEvents')
	f.create_dataset('SortInd', data=Ind, dtype=int)#to find spikes and raw data in the .txt files
	f.create_dataset('Amplitudes', data=SpkAmp[PreSelectedEvents][Ind])
	f.create_dataset('Channels', data=SpkCh[PreSelectedEvents][Ind])
	f.create_dataset('Times', data=SpkT[PreSelectedEvents][Ind])
	f.create_dataset('RepolarizingSpikes', data=PeakL[PreSelectedEvents][Ind] ,dtype=bool)
	f.create_dataset('RecalibrationOffsets', data=RecalibOffset[PreSelectedEvents][Ind])
	f.create_dataset('NumberOfCutoutChannels', data=Ntraces[PreSelectedEvents][Ind])
	g.create_dataset('NoisyChannels', data=IgnCh, dtype=int)
	f.create_dataset('PreSelectedEvents', data=PreSelectedEvents, dtype=bool)
	g.close()
	return NSpk

### Spike sorting
def readShapesFile(TxtFile, HdfFile, NSpk, Sampling, Lspike=5):
	Ncut=int(Sampling)/1002+int(Sampling)/835 +7 #length of the cutouts
	NcutL=int(Sampling)/501+int(Sampling)/835 +7 #length of the long cutouts
	Npeak=6#5#define where the peak should be (needed for subtraction of median voltage)
	#Lspike=5#max. length of the spike
	Lmax=2*int(Sampling)/1002#interval to look for sudden voltage jumps
	n=0
	f=h5py.File(HdfFile,'r+')
	g=f['RawEvents']
	g.create_dataset('Locations', (NSpk,2))
	g.create_dataset('ShAmp0', (NSpk,))
	g.create_dataset('ShAmp', (NSpk,))
	g.create_dataset('ShAmpX', (NSpk,))
	g.create_dataset('ShArea', (NSpk,))
	g.create_dataset('Shapes', (NSpk,NcutL), fillvalue=0)
	#IgnSpk=np.array(g['PreSelectedEvents'].value,dtype=bool)
	IgnCh=np.array(f['NoisyChannels'].value,dtype=int)
	SortInd=np.array(g['SortInd'].value,dtype=int)
	SpkT=np.array(g['Times'].value)
	Qdavg=np.array(f['ChannelVariability'].value)
	recCh=np.array(f['recordedChannels'].value)
	
	PeakL=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	h=f['GlobalVoltageFluctuations']
	B=3*h['medianVoltage'].value+6141#np.concatenate((3*h['medianVoltage'].value,np.zeros(20,dtype=int)+6141))
	Vsbias=h['VFBias'].value
	#invert SortInd
	SortIndInv=np.zeros(SortInd.shape[0],dtype=int)
	SortIndInv[SortInd]=np.arange(SortInd.shape[0])
	### read shapes
	#always take cumulative sum of spike shapes, after substraction of mean
	Map9=[np.array([5,1,6,2,7,3,8,4]),np.array([6,2,0,4,5]),np.array([1,6,7,3,0]),\
	np.array([4,0,2,7,8]),np.array([5,1,0,3,8]),np.array([1,0,4]),\
	np.array([2,0,1]),np.array([0,2,3]),np.array([4,0,3])]
	Map9L=np.array([8,5,5,5,5,3,3,3,3])
	Map12=[np.array([4,5,1,2,3,10,11]),np.array([4,5,6,7,2,3,0]),np.array([0,1,3,6,7,8,9]),\
	np.array([0,1,2,8,9,10,11]),np.array([0,1,11,5]),np.array([0,1,4,6]),\
	np.array([1,2,5,7]),np.array([2,1,6,8]),np.array([2,3,7,9]),\
	np.array([2,3,8,10]),np.array([0,3,11,9]),np.array([10,0,3,4])]
	Map12L=np.array([7,7,7,7,4,4,4,4,4,4,4,4])
	A9=np.array([[0,0,1,0,-1,-1,1,1,-1],[0,-1,0,1,0,-1,-1,1,1]])
	A12=np.array([[0,1,1,0,0,1,2,2,1,0,-1,-1],[0,0,1,1,-1,-1,0,1,2,2,1,0]])
	fName=np.array(['_Shapes.txt','_ShapesX.txt'])
	BInd=np.zeros((2,int(Sampling)/1002+int(Sampling)/835+2),dtype=int)#maybe a bit long...
	BInd[:,:int(Sampling)/1002-3]=1
	BInd[:,int(Sampling)/1002-3:int(Sampling)/1002]=np.arange(2,5,dtype=int)[None,:]
	BInd[0,-(int(Sampling)/835)-2:]=np.arange(NcutL-(int(Sampling)/835),NcutL+2,dtype=int)
	BInd[1,-(int(Sampling)/835)-2:]=np.arange(Ncut-(int(Sampling)/835),Ncut+2,dtype=int)
	for iiii in range(2):
		b=file(TxtFile + fName[iiii])
		for i in b:
			z=np.array(i.split(),dtype=float)
			if not ((iiii==0)*(int(z[0]) in IgnCh)):
				nInd=SortIndInv[n]
				z=1.*np.reshape(z,(-1,Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])+2))
				z[:,1]/=64.# was scaled version of 33 percentile of voltage
				#bad stuff #ignore channels in list
				goodCh=(np.min(z[:,2:]%4048,axis=1)>=48)#working channels
				badCh=np.nonzero((np.min(z[:,2:]%4048,axis=1)<48)+1*(np.in1d(z[:,0],IgnCh))+1*(z[:,0]==-1))[0]#not working or nonexistent channels
				goodChInd=np.array(z[goodCh,0],dtype=int)#indices for determining global fluctuations
				noCh=np.nonzero(z[:,0]==-1)[0]#channels that do not exist
				###global fluctuations to subtract
				bA=B[SpkT[nInd]-Npeak:SpkT[nInd]-Npeak+Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]/3.
				z[:,2:]-=np.tile(bA,(z.shape[0],1))#remove global changes
				#fine tuning based on correlations with individual channels
				bA-=np.mean(bA)
				z[goodCh,2:]-=np.outer((Vsbias[goodChInd]),bA)#remove global changes II
				z[goodCh,1]-=Vsbias[goodChInd]*bA[Npeak]#should change that in the detection instead...no! local vs. global measure...
				###normalize by variance
				z[goodCh,1:]/=np.clip(Qdavg[goodChInd][:,None],1,1e12)
				#channels (relative locations)
				zChannels=np.vstack((np.array(recCh[np.array(z[:,0],dtype=int)],dtype=int)%64,\
				np.array(recCh[np.array(z[:,0],dtype=int)],dtype=int)/64))
				if len(noCh):
					if iiii:
						zChannels[:,noCh]=A12[:,noCh]+zChannels[:,0][:,None]
					else:
						zChannels[:,noCh]=A9[:,noCh]+zChannels[:,0][:,None]
				###baseline
				z[:,1:]=(z[:,1:]-np.percentile(z[:,BInd[PeakL[nInd],:]],50,axis=1)[:,None])
				###bad channels, should treat channels that do not exist equally
				#set them to 0
				z[badCh,1:]=0
				###Amplitudes
				#minimum over Lspike consecutive frames
				Cz=np.cumsum(z[:,Npeak-Lspike/2-1:Npeak+(Lspike+9)/2],axis=1)
				Pos=np.argmin(np.sum(Cz[:,Lspike:]+Cz[:,Lspike-1:-1]-Cz[:,:-Lspike]-Cz[:,1:-Lspike+1],axis=0))
				SAmp=-(Cz[:,Pos+Lspike]+Cz[:,Pos+Lspike-1]-Cz[:,Pos]-Cz[:,Pos+1])#/(LSpike-1)#pos. values
				#look for sudden jumps
				SAmp-=np.max(np.abs(np.diff(z[:,2:Npeak-Lspike/2+Pos+1],axis=1)),axis=1)\
				+np.max(np.abs(np.diff(z[:,Npeak+(Lspike+1)/2+Pos:Lmax+1],axis=1)),axis=1)
				#interpolate (linear with 8 surrounding channels), normalize by 10 channels
				#do this interpolation after determining amplitudes only!!!
				if len(badCh):
					for jj in badCh:
						if (iiii==0):
							if jj<9:
								z[jj,1:]=np.mean(z[Map9[jj],1:],axis=0)*Map9L[jj]/8.
								SAmp[jj]=np.median(SAmp[Map9[jj]])*Map9L[jj]/8.
							else:
								z[jj,1:]=0
						else:
							if jj<12:
								z[jj,1:]=np.mean(z[Map12[jj],1:],axis=0)*Map12L[jj]/8.
								SAmp[jj]=np.median(SAmp[Map12[jj]])*Map12L[jj]/8.
							else:
								z[jj,1:]=0
				SAmp0=np.clip(SAmp,0,1e12)/2.#raw amplitude estimate
				###remove 20 percentile
				SAmp-=np.percentile(SAmp0,20)#otherwise, I might effectively add stuff...
				#clip to positive values
				SAmp=np.clip(SAmp,0,1e12)#may need to divide by 4 frames (temporal window)
				###Neighboring amplitudes
				#want to look at neighboring amplitudes and have any amplitude smaller than the sum
				#of its neighbors
				NN=np.max(np.abs(zChannels[:,:,None]-zChannels[:,None,:]),axis=0)<=1
				Oamp=np.clip(SAmp,0,np.sum(SAmp[:,None]*NN,axis=0)\
				/(1.+1.*(np.arange(z.shape[0])>=(1-iiii))+1.*(np.arange(z.shape[0])>(4-iiii))))
				###remove 50 percentile (i.e. 37.5 percentile from former baseline)
				ShAmp0=np.sum(SAmp0)
				ShAmp=np.sum(Oamp)
				if (iiii==0):
					CAmp[0]=SAmp[0]#want that noisy channels are spatially limited (but would also bias results...)
				CMamp=np.clip(Oamp-np.percentile(Oamp,50),0,1e12)
				CM=np.sum(CMamp[None,:]*zChannels,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
				#set channels with a distance larger than sqrt(2) to zero/half?
				CMamp2=np.clip(CMamp,0,np.clip(CMamp*(2-np.sqrt(np.sum((CM[:,None]-zChannels)**2,axis=0))),0,1))
				CM=np.sum(CMamp2[None,:]*zChannels,axis=1)/np.clip(np.sum(CMamp2),1e-6,1e12)
				CMamp3=np.clip(CMamp,-0.1*CMamp,CMamp*np.clip((2-np.sqrt(np.sum((CM[:,None]-zChannels)**2,axis=0))),-0.1,1))
				if np.sum(CMamp3)>np.sum(CMamp)*0.5:
					CMamp=CMamp3
					CM=np.sum(CMamp[None,:]*zChannels,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
				else:
					CMamp=np.clip(CMamp,0,CMamp*np.clip((2-np.sqrt(np.sum((CM[:,None]-zChannels)**2,axis=0))),0,1))
					CM=np.sum(CMamp[None,:]*zChannels,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
				#weighted sum of raw traces
				wA=np.sum(z[:,2:]*CMamp[:,None]*1./np.clip(np.sum(CMamp),1,1e12),axis=0)
				ShArea=np.sqrt(np.sum(Oamp[None,:]*(zChannels-CM[:,None])**2)\
				/np.clip(np.sum(Oamp),1e-6,1e12))
				g['Locations'][nInd,:]=CM[::-1]+0.5#want (y,x) format to be consistent with brainwave
				g['ShAmp0'][nInd]=ShAmp0
				g['ShAmp'][nInd]=ShAmp#new amplitudes
				g['ShAmpX'][nInd]=np.sum(CMamp)#might be less noisy... and will subtract a larger baseline for wide spikes
				g['ShArea'][nInd]=ShArea
				g['Shapes'][nInd,:Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]=wA
				n+=1
		b.close()
	g['Locations'][:,:]=np.clip(g['Locations'].value\
	+1e-2*scipy.randn(g['Locations'].value.shape[0],2),0.,63.999)#to avoid discretization effects
	f.create_dataset('lenCutouts', data=Ncut)
	f.create_dataset('lenCutoutsLong', data=NcutL)
	f.create_dataset('PeakIndex', data=Npeak)
	f.close()
	return

def IsolatedSpikes(HdfFile, IncludeLongSpikes=True, DFrames=2, MaxDist=1.):
	f=h5py.File(HdfFile,'r+')
	g=f['RawEvents']
	Sampling=f['Sampling'].value
	Loc=g['Locations'].value
	ShAmpX=g['ShAmpX'].value
	ShArea=g['ShArea'].value
	Times=g['Times'].value
	RepolarizingSpikes=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	RecalibrationOffsets=np.array(g['RecalibrationOffsets'].value)
	Shapes=np.array(g['Shapes'].value)
	Amplitudes=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/12.
	X=np.zeros((1001,2))
	Xa=np.zeros((1001))
	Xt=np.zeros((1001))
	if IncludeLongSpikes:
		Ispikes=np.ones(len(Times),dtype=bool)
		for i in range(500):
			X[i,:]=Loc[i,:]
			Xa[i]=Amplitudes[i]
			Xt[i]=int(Times[i])
		for i in range(500,len(Times)):
			X[i%1001,:]=Loc[i,:]
			Xa[i%1001]=Amplitudes[i]
			Xt[i%1001]=int(Times[i])
			j=(i-500)%500
			Ind=np.nonzero(((Xa-Xa[j])>0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind,:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Ispikes[i]=False
		for i in range(len(Times)-500,len(Times)):
			j=i%500
			Ind=np.nonzero(((Xa-Xa[j])>0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind,:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Ispikes[i]=False
	else:
		Ispikes=RepolarizingSpikes.copy()
		RepolarizingInd=np.nonzero(RepolarizingSpikes)[0]
		for i in range(500):
			k=RepolarizingInd[i]
			X[i,:]=Loc[k,:]
			Xa[i]=Amplitudes[k]
			Xt[i]=int(Times[k])
		for i in range(500,np.sum(RepolarizingSpikes)):
			k=RepolarizingInd[i]
			X[i%1001,:]=Loc[k,:]
			Xa[i%1001]=Amplitudes[k]
			Xt[i%1001]=int(Times[k])
			j=(i-500)%500
			Ind=np.nonzero(((Xa-Xa[j])>0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind,:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Ispikes[k]=False
		for i in range(np.sum(RepolarizingSpikes)-500,np.sum(RepolarizingSpikes)):
			k=RepolarizingInd[i]
			j=i%500
			Ind=np.nonzero(((Xa-Xa[j])>0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind,:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Ispikes[k]=False
	Amplitudes=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/12.
	if 'IsolatedSpikes' in g:
		del g['IsolatedSpikes']
	if 'Locations' in f:
		del f['Locations']
	if 'Amplitudes' in f:
		del f['Amplitudes']
	if 'Times' in f:
		del f['Times']
	if 'RepolarizingSpikes' in f:
		del f['RepolarizingSpikes']
	if 'RecalibrationOffsets' in f:
		del f['RecalibrationOffsets']
	if 'Shapes' in f:
		del f['Shapes']
	if 'IncludeLongSpikes' in f:
		del f['IncludeLongSpikes']
	g.create_dataset('IsolatedSpikes', data=Ispikes)
	f.create_dataset('Locations', data=Loc[Ispikes,:])
	f.create_dataset('Amplitudes', data=Amplitudes[Ispikes])
	f.create_dataset('Times', data=Times[Ispikes])
	f.create_dataset('RepolarizingSpikes', data=RepolarizingSpikes[Ispikes])
	f.create_dataset('RecalibrationOffsets', data=RecalibrationOffsets[Ispikes])
	f.create_dataset('Shapes', data=Shapes[Ispikes])
	f.create_dataset('IncludeLongSpikes', data=IncludeLongSpikes, dtype=bool)
	f.close()
	return

def CorrelationAnalysis1(HdfFile, NextN=100, NPoissonNoise=40, dClockspikes=4,\
removeRecalib=3, Nignore=50, probC=0.1, Res=3, RefRes=5):
	f=h5py.File(HdfFile,'r')
	Sampling=f['Sampling'].value
	Loc=f['Locations'].value
	tMax=f['tMax'].value
	nFrames=f['nFrames'].value
	Times=f['Times'].value
	ROffsets=f['RecalibrationOffsets'].value
	SingleSpk=np.array(f['RepolarizingSpikes'].value,dtype=bool)
	f.close()
	RInd=np.abs(ROffsets-(removeRecalib*int(Sampling))/2000)>((removeRecalib*int(Sampling))/2000-1)
	nRBins=RefRes*64
	nBins=Res*64
	Roffset=-0.5*(RefRes/Res)*(Res%2-1)+0.5*(Res-RefRes)/(Res)#+0.5*(RefRes%2-1)#in units of RefRes
	print Roffset
	NChannels=nBins**2
	
	#find local maxima---------------------------------+(0.5/RefRes)*(RefRes%2-1)
	RH=np.histogram2d(np.clip(Loc[RInd*SingleSpk,0],0,64)\
	,np.clip(Loc[RInd*SingleSpk,1],0,64)\
	,bins=(np.arange(nRBins+1)*1./RefRes,np.arange(nRBins+1)*1./RefRes))[0]
	#average over d=1/sqrt(2) channels (i.e.3x3), kernel: 1-(r/rmax)**2
	RHAvg=np.zeros(RH.shape)
	for i in range(RefRes):
		for j in range(i,RefRes):
			x=(i**2+j**2)*8./RefRes**2
			if x<1.:
				if j==0:
					RHAvg+=RH
				elif i==j:
					RHAvg[i:,i:]+=RH[:-i,:-i]*(1.-x)
					RHAvg[:-i,:-i]+=RH[i:,i:]*(1.-x)
					RHAvg[i:,:-i]+=RH[:-i,i:]*(1.-x)
					RHAvg[:-i,i:]+=RH[i:,:-i]*(1.-x)
				elif i==0:
					RHAvg[:,j:]+=RH[:,:-j]*(1.-x)
					RHAvg[j:,:]+=RH[:-j,:]*(1.-x)
					RHAvg[:,:-j]+=RH[:,j:]*(1.-x)
					RHAvg[:-j,:]+=RH[j:,:]*(1.-x)
				else:
					RHAvg[i:,j:]+=RH[:-i,:-j]*(1.-x)
					RHAvg[:-i,:-j]+=RH[i:,j:]*(1.-x)
					RHAvg[j:,i:]+=RH[:-j,:-i]*(1.-x)
					RHAvg[:-j,:-i]+=RH[j:,i:]*(1.-x)
					RHAvg[i:,:-j]+=RH[:-i,j:]*(1.-x)
					RHAvg[:-i,j:]+=RH[i:,:-j]*(1.-x)
					RHAvg[j:,:-i]+=RH[:-j,i:]*(1.-x)
					RHAvg[:-j,i:]+=RH[j:,:-i]*(1.-x)
	RHAvgMax=RHAvg.copy()
	#only keep maxima within r<0.5/sqrt(2) channels (3x3 bins)
	for i in range(RefRes-1):
		for j in range(i,RefRes-1):
			if np.sqrt(i**2+j**2)<=(RefRes/np.sqrt(8)):
				if j==0:
					print j
				elif i==j:
					RHAvgMax[i:,i:]*=RHAvg[:-i,:-i]<RHAvg[i:,i:]
					RHAvgMax[:-i,:-i]*=RHAvg[i:,i:]<RHAvg[:-i,:-i]
					RHAvgMax[i:,:-i]*=RHAvg[:-i,i:]<RHAvg[i:,:-i]
					RHAvgMax[:-i,i:]*=RHAvg[i:,:-i]<RHAvg[:-i,i:]
				elif i==0:
					RHAvgMax[:,j:]*=RHAvg[:,:-j]<RHAvg[:,j:]
					RHAvgMax[j:,:]*=RHAvg[:-j,:]<RHAvg[j:,:]
					RHAvgMax[:,:-j]*=RHAvg[:,j:]<RHAvg[:,:-j]
					RHAvgMax[:-j,:]*=RHAvg[j:,:]<RHAvg[:-j,:]
				else:
					RHAvgMax[i:,j:]*=RHAvg[:-i,:-j]<RHAvg[i:,j:]
					RHAvgMax[:-i,:-j]*=RHAvg[i:,j:]<RHAvg[:-i,:-j]
					RHAvgMax[j:,i:]*=RHAvg[:-j,:-i]<RHAvg[j:,i:]
					RHAvgMax[:-j,:-i]*=RHAvg[j:,i:]<RHAvg[:-j,:-i]
					RHAvgMax[i:,:-j]*=RHAvg[:-i,j:]<RHAvg[i:,:-j]
					RHAvgMax[:-i,j:]*=RHAvg[i:,:-j]<RHAvg[:-i,j:]
					RHAvgMax[j:,:-i]*=RHAvg[:-j,i:]<RHAvg[j:,:-i]
					RHAvgMax[:-j,i:]*=RHAvg[j:,:-i]<RHAvg[:-j,i:]
	RHAvgMax=(RHAvgMax).flatten()
	#minimum spike count (could ask for 0.01 Hz instead)
	RefInd=np.array(np.nonzero(RHAvgMax>min(2.*np.percentile(RHAvg,8),np.percentile(RHAvg,30)))[0],dtype=int)
	Freq=RHAvgMax[RefInd]
	Ncomp=len(RefInd)-Nignore
	print 'number of reference channels:', Ncomp
	print 'rate threshold:', 2*np.percentile(RHAvg,8)
	print 'rate median:', np.percentile(RHAvg,50)
	RefInd=RefInd[np.argsort(np.argsort(Freq))<Ncomp]
	#adjacent bins, large radius (21 bins)
	NN=(np.arange(-RefRes/2,RefRes/2+1)[:,None]+nRBins*np.arange(-RefRes/2,RefRes/2+1)[None,:]).flatten()\
	[((np.arange(-RefRes/2,RefRes/2+1)**2)[:,None]\
	+(np.arange(-RefRes/2,RefRes/2+1)**2)[None,:]).flatten()<=RefRes**2/4.]
	#print NN
	#small radius (9 bins)
	NNc=(np.arange(-RefRes/2,RefRes/2+1)[:,None]+nRBins*np.arange(-RefRes/2,RefRes/2+1)[None,:]).flatten()\
	[((np.arange(-RefRes/2,RefRes/2+1)**2)[:,None]\
	+(np.arange(-RefRes/2,RefRes/2+1)**2)[None,:]).flatten()<=RefRes**2/8.]
	#print NNc
	#overlapping areas should be assigned to cluster with higher spike count
	Freq=RHAvgMax[RefInd]
	RHAvgMax=np.zeros(nRBins**2,dtype=int)
	for j in np.argsort(Freq):
		for i in NN:
			RHAvgMax[RefInd[j]+i]=j+1
	#or the closer one...
	for j in np.argsort(Freq):
		for i in NNc:
			RHAvgMax[RefInd[j]+i]=j+1
	#assign local maxima to units, need to remove all assigned to 0
	# and those close to recalibration events later!
	Units=RHAvgMax[np.clip(np.array((Loc[:,1]+(0.5/RefRes)*(RefRes%2-1))*RefRes,dtype=int),0,nRBins-1)\
	+np.clip(np.array((Loc[:,0]+(0.5/RefRes)*(RefRes%2-1))*RefRes,dtype=int),0,nRBins-1)*nRBins]
	### all events
	#rasterization
	Spikes=np.clip(np.array(Res*(Loc[:,1])+0.5*(Res%2-1),dtype=int),0,nBins-1)\
	+np.clip(np.array(Res*(Loc[:,0])+0.5*(Res%2-1),dtype=int),0,nBins-1)*(nBins)
	##need this for later
	#SpikesS=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)\
	#+np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)*(nBins)
	#Histogram
	H=np.reshape(np.histogram(Spikes[RInd],bins=np.arange((nBins)**2+1))[0],(nBins,nBins))
	HCount=np.zeros(H.shape)
	for i in range(-(Res-1)/2,Res/2+1):
		i0=np.clip(i,0,nBins-1)
		i1=np.clip(nBins+i,0,nBins)
		i2=np.clip(-i,0,nBins-1)
		i3=np.clip(nBins-i,0,nBins)
		for j in range(-(Res-1)/2,Res/2+1):
			j0=np.clip(j,0,nBins-1)
			j1=np.clip(nBins+j,0,nBins)
			j2=np.clip(-j,0,nBins-1)
			j3=np.clip(nBins-j,0,nBins)
			HCount[i0:i1,j0:j1]+=H[i2:i3,j2:j3]
	HCount=HCount.flatten()
	
	if (Res%2==0):
		AddNoise=int(np.clip((tMax*NPoissonNoise)/3600/4,1,1e12))*4#4 because of doing interpolation
	else:
		AddNoise=int(np.clip((tMax*NPoissonNoise)/3600,1,1e12))
	NSpikes=HCount +Res**2*AddNoise
	MeanS=np.mean(NSpikes)
	ClockSpikes=np.arange(0,tMax*Sampling,dClockspikes,dtype=int)
	Ncs=int(np.ceil(tMax*Sampling/dClockspikes))
	#need those for later (telling where Poisson and real spikes should go for next analysis, no need for amplitudes)
	PSpikes=np.concatenate((np.zeros(Spikes.shape[0],dtype=bool), np.ones(NChannels*AddNoise,dtype=bool)))
	if (Res%2==0):#that doesn't seem to be general enough
		Map0=np.repeat(np.repeat(np.reshape(np.arange(64**2*Res**2),(64*Res,64*Res)),2,axis=1),2,axis=0)
		Map1=np.zeros((64*Res*2,64*Res*2))
		Map1[:-1,:-1]=Map0[1:,1:]
		Map1[-1,:]=0
		Map1[:,-1]=0
		Map0=Map0.flatten()
		Map1=Map1.flatten()
		PoissonNoise=np.repeat(np.arange(4*NChannels),AddNoise/4)[np.argsort(scipy.rand(AddNoise*(NChannels)))]
		Spikes=np.concatenate((Spikes,Map1[PoissonNoise],np.zeros(Ncs,dtype=int)))
		SpikesS=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)\
		+np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)*(nBins)
		SpikesS=np.concatenate((SpikesS,Map0[PoissonNoise],np.zeros(Ncs,dtype=int)))
	else:
		PoissonNoise=np.repeat(np.arange(NChannels),AddNoise)[np.argsort(scipy.rand(AddNoise*(NChannels)))]
		Spikes=np.concatenate((Spikes,PoissonNoise,np.zeros(Ncs,dtype=int)))
		#SpikesS=np.concatenate((SpikesS,PoissonNoise,np.zeros(Ncs,dtype=int)))
	Units=np.concatenate((Units,np.zeros(NChannels*AddNoise,dtype=int),np.zeros(Ncs,dtype=int)-1))
	RInd=np.concatenate((RInd,np.ones(NChannels*AddNoise,dtype=bool),np.ones(Ncs,dtype=bool)))
	SpikesT=np.concatenate((Times,\
	np.interp(scipy.rand(AddNoise*(NChannels))*(len(Times)+Ncs-1),np.arange(len(Times)+Ncs)\
	,np.sort(np.concatenate((Times,ClockSpikes))))\
	, ClockSpikes))#np.array(scipy.rand(AddNoise*(NChannels))*Sampling*tMax,dtype=int),
	SortedT=np.argsort(SpikesT)
	Spikes=Spikes[SortedT]
	Units=Units[SortedT]
	RInd=RInd[SortedT]
	SpikesT=np.sort(SpikesT)
	#remove events near recalibrations (for correlations)
	SpikesR=Spikes[RInd]
	UnitsR=Units[RInd]
	UnitsRI=np.nonzero(UnitsR)[0]#exclude Poissonspikes
	UnitsR=np.clip(UnitsR[UnitsRI],0,len(RefInd))
	#events with recalibrations for correlation index
	UnitsXI=np.nonzero(Units)[0]#exclude Poissonspikes
	UnitsX=np.clip(Units[UnitsXI],0,len(RefInd))
	#Correlation matrix
	X=np.zeros((NChannels,Ncomp+1),dtype=int)
	NZero=np.zeros((NChannels))
	k=NextN
	for i in range(UnitsRI[NextN],UnitsRI[-NextN]):
		if UnitsRI[k]<i:
			k+=1
		a=np.unique(UnitsR[k-NextN:k+NextN+1])
		b=2*NextN-len(a)+1
		if (((SpikesR[i]-(Res-1)/2)%nBins)<nBins-Res+1) and (((SpikesR[i]/nBins)-(Res-1)/2)%nBins<nBins-Res+1):
			for j in range(-(Res-1)/2,Res/2+1):
				for jj in range(-(Res-1)/2,Res/2+1):
					X[SpikesR[i]+j+jj*nBins,a]+=1
					NZero[SpikesR[i]+j+jj*nBins]+=b
		else:
			X[SpikesR[i],a]+=1
			NZero[SpikesR[i]]+=b
	print 'fraction of Clockspikes:', np.sum(NZero[1:])*1./np.sum(X)
	NSpikesR=np.sum(X,axis=0)
	NSpikesR[0]=0
	TotSpikes=np.sum(NSpikesR[1:])
	MeanSR=np.mean(NSpikesR[1:])
	print 'frequency of reference spikes:', TotSpikes*1./Sampling*1./tMax/2./NextN*1./Res**2
	print 'frequency of all spikes:', np.sum(NSpikes)*1./Sampling*1./tMax*1./Res**2
	SpikeHisto=np.zeros((NChannels,Ncomp+1))
	for i in range(1,NChannels):
		RHbin=((i/nBins)*RefRes/Res)*nRBins+((i%nBins)*RefRes/Res)
		SpikeHisto[i,:]=NSpikes[i]*NSpikesR*NextN*2./(TotSpikes-NSpikesR[RHAvgMax[RHbin]])
	#remove connections between units with different activities
	Z=np.zeros((NChannels,Ncomp+1),dtype=bool)
	for i in range(1,NChannels):#compensate too much!?
		SpikeHistoTest=(SpikeHisto[i,:]*(1.-(NZero[i]/2./NextN/(NSpikes[i]+(NSpikes[i]==0)))))
		SpikeHistoTest[SpikeHistoTest<=0.01]=0.01
		Z[i,:]=(scipy.stats.poisson.sf(X[i,:]+0.5,SpikeHistoTest)<0.999999)
	Z[0,:]=False
	Z[:,0]=False
	print np.sum(Z)/(NChannels-1.)
	SpikeHisto=np.zeros((NChannels,Ncomp+1))
	for i in range(1,NChannels):
		RHbin=((i/nBins)*RefRes/Res)*nRBins+((i%nBins)*RefRes/Res)
		fracGood=(np.sum(NSpikesR[Z[i,:]]))*1./np.sum(NSpikesR[1:])
		SpikeHisto[i,:]=NSpikes[i]*NSpikesR*NextN*2./(TotSpikes-NSpikesR[RHAvgMax[RHbin]])/fracGood
	#units with increasing # of connections
	NCorr=np.zeros((NChannels,Ncomp+1),dtype=bool)
	a=RefInd.copy()#np.arange(1,1000,dtype=int)
	for i in range(1,NChannels):
		b=np.arange(1,Ncomp+1)[(((a/nRBins*nBins-i/nBins*nRBins+Roffset*nBins)**2\
		+(a%nRBins*nBins-i%nBins*nRBins+Roffset*nBins)**2)>nBins**2*RefRes**2)]#min. distance 1 electrode
		SpikeHistoTest=(SpikeHisto[i,b]*(1.-(NZero[i]/2./NextN/(NSpikes[i]+(NSpikes[i]==0)))))
		SpikeHistoTest[SpikeHistoTest<=0.001]=0.001#ignoriere einzelne Spikes
		A=(scipy.stats.poisson.sf(X[i,b]-0.5,SpikeHistoTest)<\
		(probC*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1./np.clip(NSpikesR[b],1,1e12),1e-12,10))))
		B=(scipy.stats.poisson.sf(X[i,b]+0.5,SpikeHistoTest)<\
		(probC*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1./np.clip(NSpikesR[b],1,1e12),1e-12,10))))*(X[i,b]>0)
		Q=(scipy.stats.poisson.sf(X[i,b]-0.5,SpikeHistoTest)-\
		scipy.stats.poisson.sf(X[i,b]+0.5,SpikeHistoTest))[B-A]
		A[B-A]=scipy.rand(np.sum(B-A))<(Q/2.)
		NCorr[i,b]=A
	NCorr[0,:]=False
	NCorr[:,0]=False
	print np.sum(NCorr)/(NChannels)
	###
	#Indnc1=(np.reshape(np.arange(nBins**2),(nBins,nBins))[::2,::2]).flatten()
	#Indnc2=(np.reshape(np.arange(nBins**2),(nBins,nBins))[1::2,::2]).flatten()
	#NCorrClusterLinkage=scipy.cluster.hierarchy.ward(NCorr[Indnc1]+NCorr[Indnc1+1]\
	#+NCorr[Indnc2]+NCorr[Indnc2+1])
	#NCorrFCluster=scipy.cluster.hierarchy.fcluster(NCorrClusterLinkage,11,'maxclust')
	NCorrX=np.zeros((NChannels,Ncomp+1),dtype=bool)
	#for k in range(len(probC)):
	for i in range(1,NChannels):
		b=np.arange(1,Ncomp+1)[(((a/nRBins*nBins-i/nBins*nRBins+Roffset*nBins)**2\
		+(a%nRBins*nBins-i%nBins*nRBins+Roffset*nBins)**2)<nBins**2*nRBins**2/4)\
		*(((a/nRBins*nBins-i/nBins*nRBins+Roffset*nBins)**2\
		+(a%nRBins*nBins-i%nBins*nRBins+Roffset*nBins)**2)>nBins**2*RefRes**2)]
		SpikeHistoTest=(SpikeHisto[i,b]*(1.-(NZero[i]/2./NextN/(NSpikes[i]+(NSpikes[i]==0)))))
		SpikeHistoTest[SpikeHistoTest<=0.001]=0.001
		A=(scipy.stats.poisson.sf(X[i,b]-0.5,SpikeHistoTest)<\
		(0.25*probC*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1./np.clip(NSpikesR[b],1,1e12),1e-12,10))))
		B=(scipy.stats.poisson.sf(X[i,b]+0.5,SpikeHistoTest)<\
		(0.25*probC*np.sqrt(np.clip(MeanS*MeanSR*1./NSpikes[i]*1./np.clip(NSpikesR[b],1,1e12),1e-12,10))))*(X[i,b]>0)
		Q=(scipy.stats.poisson.sf(X[i,b]-0.5,SpikeHistoTest)-\
		scipy.stats.poisson.sf(X[i,b]+0.5,SpikeHistoTest))[B-A]
		A[B-A]=scipy.rand(np.sum(B-A))<(Q/2.)
		NCorrX[i,b]=A
	NCorrX[0,:]=False
	NCorrX[:,0]=False
	print np.sum(NCorrX)/(NChannels)
	#do a smoothening: if NCorrX==True in at least three neighbors...
	### Correlation index C
	if (Res%2==0):
		Spikes=SpikesS[SortedT]#move to regular binning
	C=np.zeros(len(Spikes))
	XY=np.concatenate((np.array([-1]),RefInd))
	k=NextN-1
	for i in range(UnitsXI[NextN],UnitsXI[-NextN]):
		if UnitsXI[k]<i:
			k+=1
			x0=np.array(UnitsX[k-NextN:k+NextN+1],dtype=int)
			y=np.array(Spikes[UnitsXI[k-NextN:k+NextN+1]],dtype=int)#mapping to pixels for Cmat
		if Spikes[i]<>0:
			a=int(Spikes[i])
			x=x0*Z[a,x0]#correlated reference channels
			a1=np.array(np.unique(x),dtype=int)[1:]
			a2=np.array(np.unique(x[NextN/2:-NextN/2+1]),dtype=int)[1:]
			a3=a1[(((a/nBins*nRBins/nBins-XY[a1]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a1]%nRBins-0.5)**2)<nRBins**2/16)\
			*(((a/nBins*nRBins/nBins-XY[a1]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a1]%nRBins-0.5)**2)>16)]
			a4=a2[(((a/nBins*nRBins/nBins-XY[a2]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a2]%nRBins-0.5)**2)<nRBins**2/16)\
			*(((a/nBins*nRBins/nBins-XY[a2]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a2]%nRBins-0.5)**2)>16)]
			a1=a1[(((a/nBins*nRBins/nBins-XY[a1]/nRBins)**2\
			+(a%nBins*nRBins/nBins-XY[a1]%nRBins-0.5)**2)<nRBins**2/2)\
			*(((a/nBins*nRBins/nBins-XY[a1]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a1]%nRBins-0.5)**2)>16)]
			a2=a2[(((a/nBins*nRBins/nBins-XY[a2]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a2]%nRBins-0.5)**2)<nRBins**2/2)\
			*(((a/nBins*nRBins/nBins-XY[a2]/nRBins-0.5)**2\
			+(a%nBins*nRBins/nBins-XY[a2]%nRBins-0.5)**2)>16)]
			a5=np.unique(y[np.in1d(x,a1)])
			a6=np.unique(y[np.in1d(x,a2)])
			a7=np.unique(y[np.in1d(x,a3)])
			a8=np.unique(y[np.in1d(x,a4)])
			c1=0
			c2=0
			c3=0
			c4=0
			if len(a4)>4:
				b=a4[NCorrX[a,a4]]
				c=a8#[NCorrX[a8,aa]]
				if len(b)>=5:
					c4=(np.sum(NCorrX[c,:][:,b])-len(c)+1)*(1./len(a4))/np.clip(len(c),1,1e6)
			if len(a3)>4:
				b=a3[NCorrX[a,a3]]
				c=a7#[NCorrX[a7,aa]]
				if len(b)>=5:
					c3=(np.sum(NCorrX[c,:][:,b])-len(c)+1)*(1./len(a3))/np.clip(len(c),1,1e6)
			if len(a2)>4:
				b=a2[NCorr[a,a2]]
				c=a6#[NCorrX[a6,aa]]
				if len(b)>=4:
					c2=(np.sum(NCorr[c,:][:,b])-len(c))*(1./len(a2))/np.clip(len(c),1,1e6)
			if len(a1)>3:
				b=a1[NCorr[a,a1]]
				c=a5#[NCorrX[a5,aa]]
				if len(b)>=4:
					c1=(np.sum(NCorr[c,:][:,b])-len(c))*(1./len(a1))/np.clip(len(c),1,1e6)
			C[i]=max(C[i],c1,c2,c3,c4)
	g=h5py.File(HdfFile,'r+')
	#SpikesAmpS=np.zeros(Spikes.shape)
	SpikesS=np.zeros(Spikes.shape,dtype=int)
	SpikesCS=np.zeros(Spikes.shape)
	SpikesTS=np.zeros(Spikes.shape)
	#RIndS=np.zeros(Spikes.shape)
	#SpikesAmpS[SortedT]=SpikesAmp#reversing sorting...
	SpikesTS[SortedT]=SpikesT
	SpikesCS[SortedT]=C
	SpikesS[SortedT]=Spikes
	#RIndS[SortedT]=RInd
	#ShInd=(SpikesAmpS[SpikesAmpS<>0]>0.05)
	#SingleSpkS=np.zeros(ShInd.shape,dtype=bool)
	#SingleSpkS[ShInd]=SingleSpk
	L=PSpikes.shape[0]
	print np.mean(C)
	print np.sum(C==0)
	print np.sum(C<>0)
	if 'CorrelationAnalysis' in g:
		del g['CorrelationAnalysis']
	if 'Units' in g:
		del g['Units']
	if 'Corr' in g:
		del g['Corr']
	f=g.create_group('CorrelationAnalysis')
	#f.create_dataset('Amplitudes', data=SpikesAmpS[SpikesAmpS<>0])
	#f.create_dataset('Cluster', data=NCorrFCluster)
	g.create_dataset('Units', data=SpikesS[:L][True-PSpikes])
	f.create_dataset('AllUnits', data=SpikesS[:L])
	f.create_dataset('PoissonUnits', data=SpikesS[:L][PSpikes])
	#f.create_dataset('Times', data=SpikesTS[:L][True-PSpikes])#no need
	f.create_dataset('PoissonTimes', data=SpikesTS[:L][PSpikes])
	g.create_dataset('Corr', data=SpikesCS[:L][True-PSpikes])
	f.create_dataset('PoissonCorr', data=SpikesCS[:L][PSpikes])
	f.create_dataset('PoissonInd', data=PSpikes)
	if 'Parameter' in f:
		del f['Parameter']
	h=f.create_group('Parameter')
	h.create_dataset('LocMax', data=RefInd)
	h.create_dataset('MaximaResolution', data=nRBins)
	h.create_dataset('Resolution', data=nBins)
	h.create_dataset('NextN', data=NextN)
	h.create_dataset('NPoissonNoise', data=NPoissonNoise)
	h.create_dataset('dClockspikes', data=dClockspikes)
	h.create_dataset('removeRecalib', data=removeRecalib)
	h.create_dataset('Nignore', data=Nignore)
	#h.create_dataset('IncludeLongSpikes', data=IncludeLongSpikes, dtype=bool)
	#f.create_dataset('RecalibFreeInd', data=RIndS[SpikesAmpS<>0])
	#f.create_dataset('ShInd', data=ShInd)#is actually redundant, but maybe saves work later... no need!
	#f.create_dataset('RepolarizingSpikes', data=SingleSpkS)#redundant, but needed in second step
	g.close()
	return

def CorrelationAnalysis2(HdfFile):
	f=h5py.File(HdfFile,'r')
	g=f['CorrelationAnalysis']
	Sampling=f['Sampling'].value
	tMax=f['tMax'].value
	nFrames=f['nFrames'].value
	CS=f['Corr'].value
	CPoisson=g['PoissonCorr'].value
	PoissonInd=g['PoissonInd'].value
	nBins=g['Parameter/Resolution'].value
	SpikesS=f['Units'].value
	SpikesPoisson=g['PoissonUnits'].value
	SpikesAll=g['AllUnits'].value
	SpikesTS=f['Times'].value
	SpikesAmpS=np.clip(f['Amplitudes'].value,1e-2,1e6)
	SingleSpk=np.array(f['RepolarizingSpikes'].value,dtype=bool)
	IncludeLongSpikes=np.array(f['IncludeLongSpikes'].value,dtype=bool)
	f.close()
	Res=nBins/64
	LIndS=np.zeros((SpikesS.shape[0],1+IncludeLongSpikes),dtype=bool)
	LInd=np.zeros((SpikesAll.shape[0],1+IncludeLongSpikes),dtype=bool)
	LIndS[:,0]=SingleSpk#short spikes
	LInd[True-PoissonInd,0]=SingleSpk#short spikes
	if IncludeLongSpikes:
		LIndS[:,1]=(True-SingleSpk)#long spikes
		LInd[True-PoissonInd,1]=(True-SingleSpk)#long spikes
	NChannels=nBins**2#np.sum(MaskX<>0)+1
	Pnew=np.zeros(len(CS.flatten()))
	NSpikes=np.histogram2d(SpikesS,(True-SingleSpk),bins=(np.arange(NChannels+1),np.arange(3)))[0]
	AmpCiKink=np.zeros((NChannels,1+IncludeLongSpikes))
	Pval=np.zeros((NChannels,1+IncludeLongSpikes))+1.
	CINoise=np.zeros((NChannels,1+IncludeLongSpikes))+1.
	fNoise=np.zeros((NChannels,1+IncludeLongSpikes))+1.
	for ii in range(1+IncludeLongSpikes):
		SameSpikes=(PoissonInd-LInd[:,ii])<>0
		for i in range(NChannels):
			if ((i%1000)==0):
				print i
			#make a list of channels
			Xind=np.clip(np.arange(i%nBins-Res/2,i%nBins+(Res+1)/2,dtype=int),0,nBins-1)
			Yind=np.clip(np.arange(i/nBins-Res/2,i/nBins+(Res+1)/2,dtype=int),0,nBins-1)
			ChListI=np.unique(Xind[:,None]+Yind[None,:]*nBins)
			if (NSpikes[i,ii])>5:
				SpikesAllInd=np.in1d(SpikesAll,ChListI)
				#SpikesAllI=(SpikesAll[SpikesAllInd]==i)
				SpikesSInd=np.in1d(SpikesS,ChListI)
				SpikesPoissonInd=np.in1d(SpikesPoisson,ChListI)
				#SpikesPoissonI=(SpikesPoisson[SpikesPoissonInd]==i)
				PoissonIndx=PoissonInd[SpikesAllInd][SameSpikes[SpikesAllInd]]
				iInd=True-PoissonIndx
				ChAmp=np.zeros(iInd.shape[0])
				ChC=np.zeros(iInd.shape[0])
				SpikesSI=np.zeros(iInd.shape[0])
				ChAmp[iInd]=SpikesAmpS[(SpikesSInd)*LIndS[:,ii]]
				ChC[iInd]=CS[(SpikesSInd)*LIndS[:,ii]]
				ChC[True-iInd]=CPoisson[(SpikesPoissonInd)]
				SpikesSI[iInd]=1.5*(SpikesS[(SpikesSInd)*(LIndS[:,ii])]==i)+0.5#weighting (surrounding spikes)
				#sort CI for Poisson spikes
				N=ChC.shape[0]
				cP=CPoisson[SpikesPoissonInd].flatten()+1e-6*scipy.rand(np.sum(SpikesPoissonInd))
				aP=np.argsort(cP)
				#sort CI for putative events
				Rnd=1e-6*scipy.rand(N)
				aX=np.argsort(ChC.flatten()+Rnd)
				aXS=np.argsort(np.argsort(aX[iInd]))
				iX=np.interp(ChC.flatten()+Rnd\
				,np.concatenate((np.array([0]),cP[aP])),np.linspace(0,1,len(aP)+1))
				iXS=iX[iInd]
				#sort amplitudes for putative events
				Rnd=1e-6*scipy.rand(N)
				aY=np.argsort(ChAmp.flatten()+Rnd)
				aYS=np.argsort(np.argsort(aY[iInd]))
				rY=np.argsort(np.argsort(ChAmp.flatten()+Rnd))
				#cumulative sums for KStest (is that one allowed?)
				X=1.*np.cumsum(rY[aX])
				X*=0.5/np.clip(X[-1],1,1e12)
				Y=np.cumsum(iX[aY])
				Y*=0.5/np.clip(Y[-1],1,1e12)
				#NS=np.sum((SpikesSInd)*(LIndS[:,ii]))
				NS=np.sum(SpikesSI)
				#Ncount[i]=NS
				#AmpCiArea[i]=np.min(X+Y-(np.arange(1,N+1)*1./N))
				#need a distortion there!
				PmapX=np.clip(np.cumsum(0.5*SpikesSI[aX]+0.5*SpikesSI[aY])-1,0,1e12)
				#SpikesS[SpikesSInd]==i
				#map for excluding Poisson spikes???
				#PmapX=np.arange(N)-np.cumsum(0.5*PoissonIndx[aX]+0.5*PoissonIndx[aY])#how many spikes on left
				BiasC=np.clip(np.sqrt(np.linspace(2./N,2,N)*np.linspace(2,2./N,N)),0.1,1)
				AmpCiKink[i,ii]=np.mean(PmapX[np.argsort((X+Y-np.arange(1,N+1)*1./N)*1./BiasC)][:NS/20])*1./(PmapX[-1])
				aiXS=np.argsort(iXS)
				CINoise[i,ii]=1.-np.mean(np.clip((iXS[aiXS]-np.cumsum(SpikesSI[iInd][aiXS])*1./NS)\
				*1./np.cumsum(SpikesSI[iInd][aiXS][::-1]+1)*NS,0,1))
				#KStest
				Pval[i,ii]=scipy.stats.kstest(X+Y, 'uniform', alternative='greater', mode='approx')[1]
				if Pval[i,ii]<0.01:
					fNoise[i,ii]=AmpCiKink[i,ii]
					#print i, Pval[i,ii]
				else:
					fNoise[i,ii]=CINoise[i,ii]
				#have to decide which spikes are noisy
				#know the fraction of noise and the distribution of noise (Poisson spikes)
				#need to average; use reflecting boundaries
				NPoisson=len(aP)
				#have (1-fNoise) space to increase Poisson probability, or as factor, mean(pPoisson)/fNoise
				#should be bound on pPoisson
				Cp=np.argsort(np.concatenate((CPoisson[SpikesPoissonInd].flatten()\
				,CS[(SpikesSInd)*LIndS[:,ii]].flatten())))<NPoisson#changed def
				Cpc=np.cumsum(np.concatenate((np.array([0]),Cp[::-1],Cp,Cp[::-1])))
				Ind=np.nonzero(Cp==0)[0]#spikes
				jj=0.025#min(0.025,0.025*int(40/np.sqrt(N)))
				w=True
				SpikesSiX=SpikesSI[iInd][aXS]
				while jj<1. and w:
					jj+=0.025
					nAll=int(N*jj)
					pPoisson=(Cpc[Ind+N+nAll+1]-Cpc[Ind+N-nAll])/(2.*nAll+1.)+1e-6
					w=(np.max(pPoisson)*fNoise[i,ii])>1.1*np.sum(pPoisson*SpikesSiX)/NS#should allow for being 10% off
				#print jj
				#
				pPoisson=(Cpc[Ind+N+nAll+1]-Cpc[Ind+N-nAll])/(2.*nAll+1.)+1e-6
				pPoisson=np.clip(pPoisson*fNoise[i,ii]*NS/np.sum(pPoisson*SpikesSiX),0,1)#does not go well with clipping...
				pPoisson=np.clip(pPoisson*fNoise[i,ii]*NS/np.sum(pPoisson*SpikesSiX),0,1)#do a round of adjustment to fNoise
				pPoisson=np.clip(pPoisson*fNoise[i,ii]*NS/np.sum(pPoisson*SpikesSiX),0,1)#do a round of adjustment to fNoise
				fTP0=np.zeros(np.sum((SpikesSInd)*(LIndS[:,ii])))
				fTP0[aXS]=1.-pPoisson
				#print np.mean(pPoisson), fNoise[i,ii]
				#now have to do the same with amplitudes
				#Amplitude_CI dependence
				#NS=np.sum((SpikesSInd)*(LIndS[:,ii]))
				'''
				nS=int(jj*np.sum((SpikesSInd)*(LIndS[:,ii])))
				AmpRefl=np.cumsum(np.concatenate((np.array([0]),fTP0[aYS][:nS][::-1],fTP0[aYS],fTP0[aYS][-nS:][::-1])))
				#normalize to one
				SpikesSiY=SpikesSI[iInd][aYS]
				pAmp=(AmpRefl[(2*nS+1):]-AmpRefl[:-(2*nS+1)])/(2.*nS+1.)+1e-6
				pAmp*=np.sum(SpikesSiY)*1./np.sum(pAmp*SpikesSiY)#does not go well with clipping... but shouldn't have that a large effect
				#print np.mean(pAmp)
				#modify probabilities multiplicatively and locally
				Cpc2=np.cumsum(np.concatenate((np.array([0]),fTP0[:nS][::-1],fTP0,fTP0[-nS:][::-1])))
				fTP1=np.zeros(np.sum((SpikesSInd)*(LIndS[:,ii])))
				fTP1[aYS]=np.clip(fTP0[aYS]*pAmp,0,1)#will correct for clipping in the normalization step
				Cpc3=np.cumsum(np.concatenate((np.array([0]),fTP1[:nS][::-1],fTP1,fTP1[-nS:][::-1])))
				Cnew=(Cpc3[2*nS+1:]-Cpc3[:-2*nS-1])*0.5/(nS+0.5)
				Cold=(Cpc2[2*nS+1:]-Cpc2[:-2*nS-1])*0.5/(nS+0.5)
				fTP1*=Cold*1./np.clip(Cnew,1e-12,1e12)#correct for a local amplitude bias (in CI)
				#print np.mean(fTP1), np.mean(np.clip(fTP1[SpikesSI[iInd]>1],0,1))
				Pnew[(SpikesS==i)*LIndS[:,ii]]+=np.clip(fTP1[SpikesSI[iInd]>1],0,1)
				'''
				Pnew[(SpikesS==i)*LIndS[:,ii]]+=np.clip(fTP0[SpikesSI[iInd]>1],0,1)
	#save results
	g=h5py.File(HdfFile,'r+')
	f=g['CorrelationAnalysis']
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

def Clustering(HdfFile, fNextMin=0.02,  fNext=0.1, Resolution=12,\
gradientThr=0.5, MaxDist=1., pOffset=0.2, useProbabilities=True):
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	rCh=g['recordedChannels'].value
	NCh=rCh.shape[0]
	Loc=g['Locations'].value
	f=g['CorrelationAnalysis']
	Probability=f['Probability'].value
	#ISpikes=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	g.close()
	if not useProbabilities:
		pOffset=0
	nNextMin=int(fNextMin*tMax)
	nNext=int(fNext*tMax)
	nBins=Resolution*64
	#make 2d gradient map, recursively find maxima
	Dr=np.zeros((nBins**2))
	Dr100=np.zeros((nBins**2))
	Da=np.zeros((nBins**2),dtype=int)
	Qinc=np.array([0,-nBins-1,-nBins,-nBins+1,1,+nBins+1,+nBins,+nBins-1,-1])
	#this is just to have smaller arrays/matrices (split problem in parts)
	for i in range(16):
		print i
		LocY=Loc[np.abs(Loc[:,0]-4*i-2)<4,:]
		if useProbabilities:
			PY=Probability[np.abs(Loc[:,0]-4*i-2)<4]
		for j in range(16):
			LocX=LocY[np.abs(LocY[:,1]-4*j-2)<4,:]
			if useProbabilities:
				PX=np.clip(PY[np.abs(LocY[:,1]-4*j-2)<4]-pOffset,0,1.)
			else:
				PX=np.ones(LocX.shape[0])
			for k in range(nBins/16):
				for l in range(nBins/16):
					X=LocX-np.array([i*4+k*64./nBins,j*4+l*64./nBins])[None,:]
					#Distances
					Y=np.sum(X**2,axis=1)
					Yind=(Y<MaxDist)
					#sorting according to distance (for distances of less than 42 mum)
					Z=np.argsort(Y[Yind])[:nNext]
					#have to find maximum over a range of neighbors, ignore closest 50 spikes
					if len(Z)>nNextMin:
						Xnorm=X[Yind,:][Z,:]#surrounding spikes
						Pnorm=PX[Yind][Z]*np.clip(Y[Yind][Z]*Resolution**2,0,1)
						#summing unit vectors, weighted
						Q=np.cumsum(np.cumsum(Xnorm*Pnorm[:,None]*1./np.clip(np.sqrt(np.sum(Xnorm**2,axis=1)),1e-6,1e6)[:,None]\
						,axis=0),axis=0)[nNextMin:,:]\
						/np.sqrt(np.cumsum(np.cumsum(Pnorm+pOffset)))[nNextMin:,None]
						#how many spikes to look at...
						QzInd=np.argmax(np.sum(Q**2,axis=1))
						#very local behaviour (spatial spread)
						Dr100[(i*nBins/16+k)*nBins+j*nBins/16+l]=np.sum(Xnorm[nNextMin,:]**2)#why -30? removed.
						#circular variance
						Qz=np.sqrt(np.sum(Q[QzInd,:]**2))
						if Qz>gradientThr:#forgot an else here?
							Dr[(i*nBins/16+k)*nBins+j*nBins/16+l]=Qz
							Da[(i*nBins/16+k)*nBins+j*nBins/16+l]=int((np.arctan2(Q[QzInd,0],Q[QzInd,1])/np.pi*4.+3.5)%8+1)
					else:
						Dr[(i*nBins/16+k)*nBins+j*nBins/16+l]=-1
	#Clustering, start with grid
	CMap=np.arange(nBins**2,dtype=int)
	CMap[Dr==-1]=0#ignore areas with low counts (0 is garbage cluster)
	#need to cumulate clusters where there is no gradient
	Step=np.array([1,nBins,nBins+1])
	for i in range(nBins*(nBins-1)-1):
		if (i%nBins)<>(nBins-1):
			for k in range(3):
				if Da[i]==0 and Da[i+Step[k]]==0:
					CMap[i+Step[k]]=min(CMap[i],CMap[i+Step[k]])
	#do that again as I only propagate in two directions
	for i in range(nBins*(nBins-1)-1):
		if (i%nBins)<>(nBins-1):
			for k in range(3):
				if Da[i]==0 and Da[i+Step[k]]==0:
					if CMap[i]<>CMap[i+Step[k]]:
						if CMap[i]<CMap[i+Step[k]]:
							CMap[CMap==CMap[i+Step[k]]]=CMap[i]
						else:
							CMap[CMap==CMap[i]]=CMap[i+Step[k]]
	#need one step where I reassign spikes to center locations of no-gradient areas here
	for i in np.nonzero(np.bincount(CMap)>2)[0][1:]:
		q=np.nonzero(CMap==i)[0]
		NewLoc=int(np.median(q%nBins))+int(np.median(q/nBins))*nBins
		if i<>NewLoc:
			if NewLoc in q:
				CMap[q]=NewLoc
			else:
				CMap[q]=q[np.argmin((q%nBins-NewLoc%nBins)**2+(q/nBins-NewLoc/nBins)**2)]
	#need to find where directions are consistent with neighbor
	a=Qinc[Da[(Qinc[Da]+CMap)%(nBins**2)]]
	b=Qinc[Da]
	c=(a>0)*(b>0)
	d=np.abs(np.abs((a[c]-1)-(b[c]-1))-4)<3#take those only if CV greater than opposite direction
	Ind=np.nonzero(c[d])[0][Dr[c[d]]<=Dr[Qinc[Da[c[d]]]+CMap[c[d]]]+gradientThr]
	Da[Ind]=0
	for i in range(nBins*(nBins-1)-1):
		if (i%nBins)<>(nBins-1):
			for k in range(3):
				if Da[i]==0 and Da[i+Step[k]]==0:
					if CMap[i]<>CMap[i+Step[k]]:
						if CMap[i]<CMap[i+Step[k]]:
							CMap[CMap==CMap[i+Step[k]]]=CMap[i]
						else:
							CMap[CMap==CMap[i]]=CMap[i+Step[k]]
	for i in range(nBins*(nBins-1)-1):
		if (i%nBins)<>(nBins-1):
			for k in range(3):
				if Da[i]==0 and Da[i+Step[k]]==0:
					if CMap[i]<>CMap[i+Step[k]]:
						if CMap[i]<CMap[i+Step[k]]:
							CMap[CMap==CMap[i+Step[k]]]=CMap[i]
						else:
							CMap[CMap==CMap[i]]=CMap[i+Step[k]]
	#gradient descent
	Qnext=np.arange(nBins**2,dtype=int)+Qinc[Da]
	ClusterM=CMap[Qnext]
	i=0
	while (ClusterM-ClusterM[Qnext]).any():
		Ind=scipy.rand(nBins**2)>0.1
		ClusterM[Ind]=ClusterM[Qnext][Ind]
		i+=1
		if i>500:
			print 'something wrong'
			break
	print i
	#check for nearby peaks (all surrounding area belongs to another peak?)
	for i in np.unique(ClusterM):
		a=np.clip(np.array([i%nBins-1,i%nBins+2,i/nBins-1,i/nBins+2],dtype=int),0,nBins)
		b=(np.arange(a[0],a[1],dtype=int)[:,None]+nBins*np.arange(a[2],a[3],dtype=int)[None,:]).flatten()
		if ((np.sum(ClusterM[b]==ClusterM[i])-1)*1./(b.shape[0]-1))<0.5:
			c=np.argmax(np.diff(np.concatenate((np.array([-1])\
			,np.nonzero(np.diff(np.sort(ClusterM[b])))[0],np.array([b.shape[0]])))))
			ClusterM[ClusterM==i]=np.sort(np.unique(ClusterM[b]))[c]
	if nBins>=768:
		for i in np.unique(ClusterM):
			a=np.clip(np.array([i%nBins-2,i%nBins+3,i/nBins-2,i/nBins+3],dtype=int),0,nBins)
			b=(np.arange(a[0],a[1],dtype=int)[:,None]+nBins*np.arange(a[2],a[3],dtype=int)[None,:]).flatten()
			if ((np.sum(ClusterM[b]==ClusterM[i])-1)*1./(b.shape[0]-1))<0.5:
				c=np.argmax(np.diff(np.concatenate((np.array([-1])\
				,np.nonzero(np.diff(np.sort(ClusterM[b])))[0],np.array([b.shape[0]])))))
				ClusterM[ClusterM==i]=np.sort(np.unique(ClusterM[b]))[c]
	for i in np.unique(ClusterM):
		j=int(nBins>=768)
		a=np.clip(np.array([i%nBins-2-j,i%nBins+3+j,i/nBins-2-j,i/nBins+3+j],dtype=int),0,nBins)
		b=(np.arange(a[0],a[1],dtype=int)[:,None]+nBins*np.arange(a[2],a[3],dtype=int)[None,:]).flatten()
		if ((np.sum(ClusterM[b]==ClusterM[i])-1)*1./(b.shape[0]-1))<0.7:
			b=b[ClusterM[b]<>i]
			c=np.argmax(np.diff(np.concatenate((np.array([-1])\
			,np.nonzero(np.diff(np.sort(ClusterM[b])))[0],np.array([b.shape[0]])))))
			ClusterM[ClusterM==i]=np.sort(np.unique(ClusterM[b]))[c]
	ClusterList=np.sort(np.unique(ClusterM))
	CAreaMatrix=np.zeros((nBins*nBins))
	n=0
	for i in ClusterList:
		CAreaMatrix[ClusterM==i]=n
		n+=1			
	CAreaMatrix=np.array(CAreaMatrix,dtype=int)
	#new spiketrains
	xVal=np.array(Loc[:,1]*Resolution,dtype=int)
	yVal=np.array(Loc[:,0]*Resolution,dtype=int)
	SpkLocId=np.reshape(CAreaMatrix,(nBins,nBins))[yVal,xVal]
	#Count=np.histogram(SpkLocId,bins=np.arange(n+1))[0]
	#reassign clusters
	for i in range(n):
		if np.sum(SpkLocId==i)==0:#no spikes
			CAreaMatrix[CAreaMatrix==i]=0
			ClusterList[i]=0
	#remove outliers
	for i in range(1,n):
		Ind=np.nonzero(CAreaMatrix==i)[0]
		if len(Ind):
			M=np.min(Dr100[Ind])
			q=Ind[Dr100[Ind]>36*M]
			CAreaMatrix[q]=0#remove areas where density is 36 times lower than at max.
	CBoundaries=np.zeros(nBins*nBins,dtype=int)
	Dr100min=np.zeros(n)	
	Dr100peak=np.zeros(n)	
	#find shells
	for i in range(1,n):
		Ind=np.nonzero(CAreaMatrix==i)[0]
		if len(Ind):
			#outer shell
			A=np.zeros(len(Ind),dtype=int)
			for j in range(len(Ind)):
				if np.sum(np.in1d(Ind[j]+Qinc,Ind))>7:#bin is not in shell
					A[j]=1
			#next layer
			for k in range(1,3):
				Indin=np.nonzero(A==k)[0]
				for j in Indin:
					if np.sum(np.in1d(Ind[j]+Qinc,Ind[Indin]))>7:#bin is not in shell
						A[j]=k+1
			#going back, finding the adjacent locations
			for k in range(3)[::-1]:
				Indin=np.nonzero(A==k)[0]
				Indout=np.nonzero(A>k)[0]
				for j in Indin:
					if np.sum(np.in1d(Ind[j]+Qinc,Ind[Indout]))<2:#bin is not in shell
						A[j]=k-1
			#discarding bins
			CAreaMatrix[Ind[A==-1]]=0
			CBoundaries[Ind]=A
			#computing min Dr100 at boundary
			if np.sum(A==0):
				Dr100min[i]=np.min(Dr100[Ind[A==0]])
				Dr100peak[i]=np.min(Dr100[Ind[A>0]])
	SpkLocId=np.reshape(CAreaMatrix,(nBins,nBins))[yVal,xVal]
	Mapv=np.bincount(np.unique(SpkLocId))
	MapNew=np.cumsum(Mapv)-1
	nn=MapNew[-1]+1
	print nn
	for i in range(n):
		if np.sum(SpkLocId==i)==0:#no spikes
			CAreaMatrix[CAreaMatrix==i]=0
	CAreaMatrix=MapNew[CAreaMatrix]
	SpkLocId=MapNew[SpkLocId]
	ClusterList=ClusterList[np.array(Mapv,dtype=bool)]
	Dr100min=Dr100min[np.array(Mapv,dtype=bool)]
	Dr100peak=Dr100peak[np.array(Mapv,dtype=bool)]
	Count=np.histogram(SpkLocId,bins=np.arange(nn+1))[0]
	CountLoc=np.reshape(Count[CAreaMatrix],(nBins,nBins))
	CountNew=np.histogram(SpkLocId,bins=np.arange(nn+1))[0]
	Areas=np.histogram(np.histogram(CAreaMatrix,bins=np.arange(nn+1))[0],bins=np.arange(5*(nBins/128.)**2+1)*4.)[0]
	Noise=np.histogram(SpkLocId,bins=np.arange(nn+1),weights=Probability)[0]*1./np.clip(CountNew,1,1e12)
	#need to determine new variances and plot
	## I want to reduce the contribution of the noise,
	## which gives a uniform offset and increases the variance. On the 
	## other hand, I want to not bias small spike counts towards having a lower variance
	##fit a parabola and estimate the curvature? -- would only use the peak region
	##two parameters (a-bx**2) and also exact location of the peak necessary
	##variance of counts near the center? -- count small Dr100
	##compare next-neighbor-difference variance with global variability? (ratio)
	##maximize?
	##use sqrt(fourth moment) for a subtractive normalization (interested in an intermediate scale)
	V=np.zeros(nn)
	for i in range(1,nn):
		if np.sum(SpkLocId==i)>2:
			V[i]=np.sum(np.var(Loc[SpkLocId==i,:],axis=0))
			#A=np.cumsum(np.sort(CountNew[SpkLocId==i])[::-1])
			#V[i]=np.interp(0.5,A*1./A[-1],np.arange(len(A)))
	VLoc=np.reshape(V[CAreaMatrix],(nBins,nBins))
	g=h5py.File(HdfFile,'r+')
	if 'Cluster' in g:
		del g['Cluster']
	f=g.create_group('Cluster')
	f.create_dataset('ClusterId', data=SpkLocId)
	f.create_dataset('NCount', data=CountNew)
	f.create_dataset('ClusterVariance', data=VLoc)
	f.create_dataset('CLocations', data=ClusterList)
	f.create_dataset('gradientMap', data=Da)
	f.create_dataset('CAreaMatrix', data=CAreaMatrix)
	f.create_dataset('CDensityDecay', data=Dr100peak*1./np.clip(Dr100min,1e-6,1e6))
	f.create_dataset('CBoundaries', data=CBoundaries)
	f.create_dataset('Noise', data=Noise)
	h=f.create_group('Parameter')
	h.create_dataset('nBins', data=nBins)
	h.create_dataset('fNext', data=fNext)
	h.create_dataset('fNextMin', data=fNextMin)
	h.create_dataset('gradientThreshold', data=gradientThr)
	h.create_dataset('pOffset', data=pOffset)
	h.create_dataset('useProbabilities', data=useProbabilities, dtype=bool)
	g.close()
	return

def FindNoisyChannels(HdfFileList, NoisyChFile, nThreshold=3):
	#rather want to use the fraction of noise field
	#further want a more intrinsic threshold (take the average noise level and compute a factor)
	#remove everything more than 10x average noise level, but need to compensate for highly active channels
	#NoisyCh=np.zeros(4096,dtype=int)
	N=np.zeros(4096,dtype=int)
	Sp=np.zeros(4096)
	nFiles=len(HdfFileList)
	for i in range(nFiles):
		g=h5py.File(HdfFileList[i],'r')
		if i==0:
			recCh=g['recordedChannels'].value
			nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
			Sp0=np.zeros((nBins0,nBins0))
			N0=np.zeros(nBins0**2,dtype=int)
		else:
			assert len(recCh)==len(g['recordedChannels'].value)
			assert nBins0==g['CorrelationAnalysis/Parameter/Resolution'].value
		tMax=g['tMax'].value
		Loc=g['Locations'].value
		fNoise=g['CorrelationAnalysis/Noise'].value
		#P=g['CorrelationAnalysis/Probability'].value
		g.close()
		N=np.clip(np.histogram2d(Loc[:,0],Loc[:,1],bins=(np.arange(65),np.arange(65)))[0],1,1e12).flatten()
		N0=np.clip(np.histogram2d(Loc[:,0],Loc[:,1]\
		,bins=(np.arange(nBins0+1)*64./nBins0,np.arange(nBins0+1)*64./nBins0))[0],1,1e12).flatten()
		Sp0+=(np.histogram2d(Loc[:,0],Loc[:,1]\
		,bins=(np.arange(nBins0+1)*64./nBins0,np.arange(nBins0+1)*64./nBins0))[0])\
		*np.sum(np.reshape(fNoise, (nBins0,nBins0,2)),axis=2)
	for i0 in range(nBins0/64):
		for i1 in range(nBins0/64):
			Sp+=(Sp0[i0::nBins0/64,i1::nBins0/64]).flatten()
	#Pavg=Sp*1./N
	#Pavg0=Sp0*1./N0
	pThreshold=nThreshold*np.median(Sp)
	pThreshold0=nThreshold*np.median(Sp0)
	X=Sp*np.sqrt(np.median(N)*1./N)>pThreshold
	X0=Sp0.flatten()*np.sqrt(np.median(N0)*1./N0)>pThreshold0-1
	X1=np.concatenate((np.zeros((2,nBins0))\
	,np.cumsum(np.reshape(X0,(nBins0,nBins0)),axis=0),np.zeros((1,nBins0))),axis=0)
	X2=np.concatenate((np.zeros((nBins0,2))\
	,np.cumsum(X1[3:,:]-X1[:-3,:],axis=1),np.zeros((nBins0,1))),axis=1)
	X3=(X2[:,3:]-X2[:,:-3])>2#remove if less than 2 neighbors
	X4=(Sp0.flatten()*np.sqrt(np.median(N0)*1./N0)>pThreshold0)*X3.flatten()
	Y=np.arange(len(recCh))[X[recCh]]
	print Y
	print len(Y)
	print np.sum(X4)
	f=h5py.File(NoisyChFile,'w')
	f.create_dataset('NoisyChannels', data=Y)
	f.create_dataset('NoisyAreas', data=X4)
	f.close()
	return

def PeriEventActivity(HdfFile,NoisyChFile, minAmp=2,Res=12,Ns0=5./12,Na0=4, nNext=200):
	f=h5py.File(NoisyChFile,'r+')
	NoisyAreas=f['NoisyAreas'].value
	f.close()
	nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	Loc=g['Locations'].value
	Amp=g['Amplitudes'].value
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res0=int(nBins0/64)
	Ns=int(Ns0*Res)
	Na=int(Na0*Res)
	nBins=64*Res
	Spikes=np.clip(np.array(Res0*(Loc[:,1]),dtype=int),0,nBins0-1)\
	+np.clip(np.array(Res0*(Loc[:,0]),dtype=int),0,nBins0-1)*(nBins0)
	#H=np.histogram(Spikes,bins=np.arange(nBins0**2+1))[0]
	GoodCh=np.nonzero(True-NoisyAreas)[0]#np.nonzero(1-(NNoise>0.7)*(H>(tMax/Res0**2)))[0]#remove only highly noisy channels
	Mask=np.array(np.zeros(nBins0**2)-1,dtype=int)
	NChannels=len(GoodCh)
	print NChannels
	Mask[GoodCh]=np.arange(NChannels)
	#take good channels
	Loc=Loc[(Mask[Spikes]<>-1)*(Amp>minAmp),:]
	#Amp=Amp[(Mask[Spikes]<>-1),:]
	Spikes=Mask[Spikes[(Mask[Spikes]<>-1)*(Amp>minAmp)]]
	SpikesC=Loc[:,1]+1j*Loc[:,0]
	N=Spikes.shape[0]
	#make smoothened histograms of surrounding spikes
	Dhist=np.zeros((N,5))
	Ahist=np.zeros((N,6))
	for i in range(nNext,N-nNext):
		D=np.abs(SpikesC[i-nNext:i+nNext+1]-SpikesC[i])
		Dx=(D<26)*(D>2)
		Dl1=np.cumsum(Dx)
		Dl1-=Dl1[nNext]
		Ind=np.nonzero(np.abs(Dl1<31)*Dx)[0]
		Dc=np.clip(((D[Ind]-2.)/6.),0,4)
		Nc=len(Dc)
		Dhist[i,:]+=np.histogram(Dc,bins=np.arange(6),weights=1.-Dc%1)[0]*1./(Nc+(Nc==0))
		Dhist[i,:]+=np.histogram(Dc+1,bins=np.arange(6),weights=Dc%1)[0]*1./(Nc+(Nc==0))
		Ac=np.angle(SpikesC[i-nNext:i+nNext+1][Ind]-SpikesC[i])*3./np.pi-0.5
		Ahist[i,:]+=np.histogram(Ac%6,bins=np.arange(7),weights=1.-Ac%1)[0]*1./(Nc+(Nc==0))
		Ahist[i,:]+=np.histogram((Ac+1)%6,bins=np.arange(7),weights=Ac%1)[0]*1./(Nc+(Nc==0))
	#Average locally
	Dloc=np.zeros((nBins,nBins,5))
	Aloc=np.zeros((nBins,nBins,6))
	for i in range(5):
		Dloc[:,:,i]=np.histogram2d(Loc[nNext:-nNext,0]*nBins/64.,Loc[nNext:-nNext,1]*nBins/64.,\
		bins=(np.arange(nBins+1),np.arange(nBins+1)), weights=Dhist[nNext:-nNext,i])[0]
	for i in range(6):
		Aloc[:,:,i]=np.histogram2d(Loc[nNext:-nNext,0]*nBins/64.,Loc[nNext:-nNext,1]*nBins/64.,\
		bins=(np.arange(nBins+1),np.arange(nBins+1)), weights=Ahist[nNext:-nNext,i])[0]
	Nloc=np.histogram2d(Loc[nNext:-nNext,0]*nBins/64.,Loc[nNext:-nNext,1]*nBins/64.,\
	bins=(np.arange(nBins+1),np.arange(nBins+1)))[0]
	#smoothening
	H0=np.concatenate((np.zeros((1,nBins,5)),np.cumsum(Dloc,axis=0)))
	H0=H0[Ns:,:,:]-H0[:-Ns,:,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1,5)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:,:]-H1[:,:-Ns,:]
	Dloc1=H1*1./Ns**2
	H0=np.concatenate((np.zeros((1,nBins,6)),np.cumsum(Aloc,axis=0)))
	H0=H0[Ns:,:,:]-H0[:-Ns,:,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1,6)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:,:]-H1[:,:-Ns,:]
	Aloc1=H1*1./Ns**2
	H0=np.concatenate((np.zeros((1,nBins)),np.cumsum(Nloc,axis=0)))
	H0=H0[Ns:,:]-H0[:-Ns,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:]-H1[:,:-Ns]
	Nloc1=H1*1./Ns**2
	#Normalization
	Dloc1=Dloc1*1./np.clip(Nloc1,1e-6,1e12)[:,:,None]
	Aloc1=Aloc1*1./np.clip(Nloc1,1e-6,1e12)[:,:,None]
	#local averages
	X=np.concatenate((np.arange(Na)[::-1],np.arange(nBins-Ns+1),np.arange(nBins-Na-Ns+1,nBins-Ns+1)[::-1]))
	H0=np.cumsum((np.sum(Dloc1[X,:,:],axis=2)>0)*1.,axis=0)
	H0=H0[2*Na:,:]-H0[:-2*Na,:]
	H1=np.cumsum(H0[:,X],axis=1)
	Navg=np.clip(H1[:,2*Na:]-H1[:,:-2*Na],1,1e12)
	H0=np.cumsum(Dloc1[X,:,:],axis=0)
	H0=H0[2*Na:,:,:]-H0[:-2*Na,:,:]
	H1=np.cumsum(H0[:,X,:],axis=1)
	H1=H1[:,2*Na:,:]-H1[:,:-2*Na,:]
	Davg=H1*1./Navg[:,:,None]
	H0=np.cumsum(Aloc1[X,:,:],axis=0)
	H0=H0[2*Na:,:,:]-H0[:-2*Na,:,:]
	H1=np.cumsum(H0[:,X,:],axis=1)
	H1=H1[:,2*Na:,:]-H1[:,:-2*Na,:]
	Aavg=H1*1./Navg[:,:,None]
	#need to save Aavg, Aloc1, Nloc1, Davg, Dloc1 and some parameters
	if 'PeriEventActivity' in g:
		del g['PeriEventActivity']
	f=g.create_group('PeriEventActivity')
	f.create_dataset('Aavg', data=Aavg)
	f.create_dataset('Aloc', data=Aloc1)
	f.create_dataset('Davg', data=Davg)
	f.create_dataset('Dloc', data=Dloc1)
	f.create_dataset('Nloc', data=Nloc1)
	f.create_dataset('Na', data=Na)
	f.create_dataset('Ns', data=Ns)
	f.create_dataset('Res', data=Res)
	f.create_dataset('nNext', data=nNext)
	g.close()
	return
