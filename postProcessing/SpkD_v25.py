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
	#n=0
	n1=0
	b=file(TxtFile + '_Info.txt')
	for i in b:
		if '#' in i:
			k+=1
			print i
			if k==9:
				XCh4=np.zeros((len(recCh),12),dtype=int)
				n0=0
			if k==10:
				XCh5=np.zeros((len(recCh),9),dtype=int)
				SIprod=np.zeros((len(recCh),13),dtype=long)
				QdAvg=np.zeros((len(recCh)),dtype=int)
				QdAvgN=np.zeros((len(recCh)),dtype=int)
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
				cutoutDim=np.array(i.split(),dtype=int)
			if k==8:
				recCh.append(int(i))
			if k==9:
				XCh4[n0,:]=np.array(i.split(),dtype=int)
				n0+=1
			if k==10:
				XCh5[n0,:]=np.array(i.split(),dtype=int)
				n0+=1
			if k==11:
				z=np.array(i.split(),dtype=int)
				X.append(z[0])
				QdAvg+=z[1:]
				QdAvgN+=z[1:]>0
				#n+=1
			if k==12:
				SqIglobal=np.array(i.split(),dtype=long)
			if k==13:
				SqIv=np.array(i.split(),dtype=long)
			if k==14:
				SIprod[:,n1]=np.array(i.split(),dtype=long)
				n1+=1
			if k==15:
				Vsbias=-np.array(i.split(),dtype=long)*1./Sampling/nFrames#in ADCcounts
			if k==16:
				Vsqbias=np.array(i.split(),dtype=long)*1./Sampling**2/nFrames
	b.close()
	Recalib=np.array(X)# was used for recalibration events
	Qdavg=np.clip(QdAvg,640,1e12)*1./np.clip(QdAvgN,1,1e12)/64.#/n#in ADCcounts (for individual channels!)
	recCh=np.array(recCh,dtype=int)
	NCh=len(recCh)
	f=h5py.File(HdfFile,'w')
	#Parameters of the recording and spike detection
	f.create_dataset('Sampling', data=Sampling)
	f.create_dataset('reverse Detection', data=(cutoutDim[-1]<0))
	f.create_dataset('PreCut', data=cutoutDim[0])
	f.create_dataset('PostCut', data=cutoutDim[1])
	f.create_dataset('NCut', data=cutoutDim[2])
	f.create_dataset('NCutLong', data=cutoutDim[3])
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
	f.create_dataset('RecalibrationEvents', data=Recalib[::np.sign(cutoutDim[-1])])# timestamps of recalibration events
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
def readAvgFile(TxtFile, HdfFile, NCh):
	b=file(TxtFile + '_Avg.txt')
	X=[]
	for i in b:
		X.append(int(i))
	b.close()
	B=np.array(X)
	Bstd=np.std(B)/3.
	f=h5py.File(HdfFile,'r+')
	nFrames=f['nFrames'].value
	if B.shape[0]<nFrames:####may be obsolete
		print nFrames, B.shape[0]
		if 'reverse Detection' in f:
			if f['reverse Detection'].value:
				B=np.concatenate((np.ones(nFrames-B.shape[0])*B[-1],B[::-1]))
			else:
				B=np.concatenate((B,np.ones(nFrames-B.shape[0])*B[-1]))
		else:
			B=np.concatenate((B,np.ones(nFrames-B.shape[0])*B[-1]))
	g=f['GlobalVoltageFluctuations']
	g.create_dataset('medianVoltage', data=B/3.)
	g.create_dataset('medianVoltageStd', data=Bstd)
	f.close()
	#B=np.concatenate((B,np.ones(100)*B[-1]))
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
		Amp.append((z[2])*1./np.clip(z[3],1,1e12))
		Y1.append(z[4])
		YR.append(z[5])
		Y2.append(9+(np.sum((z[6:9]+z[7:])/2)+(z[6]+z[9])/2)*3)
	b.close()
	b=file(TxtFile + '_SpikesX.txt')
	for i in b:
		z=np.array(i.split(),dtype=int)
		X.append(z[0]+NCh)
		Y.append(z[1])
		Amp.append((z[2])*1./np.clip(z[3],1,1e12))
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
def readShapesFile(TxtFile, HdfFile, NSpk, Sampling=7000, Lspike=4):
	n=0
	f=h5py.File(HdfFile,'r+')
	Sampling=int(f['Sampling'].value)
	Ncut=int(f['NCut'].value)
	NcutL=int(f['NCutLong'].value)
	PreCut=int(f['PreCut'].value)
	PostCut=int(f['PostCut'].value)
	Reverse=int(f['reverse Detection'].value)
	if Reverse:
		Bsupport=(2*PreCut)/3+PostCut/3+Ncut/6
		BInd=np.zeros((2,Bsupport),dtype=int)#maybe a bit long...
		BInd[:,:Ncut/6]=1#online baseline estimate
		BInd[:,Ncut/6:Ncut/6+(2*PreCut)/3]=2+np.arange((2*PreCut)/3,dtype=int)[None,:]
		BInd[0,-PostCut/3:]=NcutL+2+np.arange(-PostCut/3,0,dtype=int)
		BInd[1,-PostCut/3:]=Ncut+2+np.arange(-PostCut/3,0,dtype=int)
	else:
		Bsupport=PreCut/3+(2*PostCut)/3+Ncut/6
		BInd=np.zeros((2,Bsupport),dtype=int)#maybe a bit long...
		BInd[:,:Ncut/6]=1#online baseline estimate
		BInd[:,Ncut/6:Ncut/6+PreCut/3]=2+np.arange(PreCut/3,dtype=int)[None,:]
		BInd[0,-(2*PostCut)/3:]=NcutL+2+np.arange(-(2*PostCut)/3,0,dtype=int)
		BInd[1,-(2*PostCut)/3:]=Ncut+2+np.arange(-(2*PostCut)/3,0,dtype=int)
	g=f['RawEvents']
	g.create_dataset('Locations', (NSpk,2))
	g.create_dataset('ShAmp0', (NSpk,))
	g.create_dataset('ShAmp', (NSpk,))
	g.create_dataset('ShAmpX', (NSpk,))
	g.create_dataset('ShArea', (NSpk,))
	g.create_dataset('Shapes', (NSpk,NcutL), fillvalue=0)
	IgnCh=np.array(f['NoisyChannels'].value,dtype=int)
	SortInd=np.array(g['SortInd'].value,dtype=int)
	SpkT=np.array(g['Times'].value)
	Qdavg=np.array(f['ChannelVariability'].value)
	recCh=np.array(f['recordedChannels'].value)
	
	PeakL=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	h=f['GlobalVoltageFluctuations']
	#B=np.concatenate((3*h['medianVoltage'].value,np.zeros(NcutL,dtype=int)+6141))
	B=3*h['medianVoltage'].value
	print B.shape
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
				bA=B[SpkT[nInd]-PreCut:SpkT[nInd]-PreCut+Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]/3.
				#print bA.shape
				bA2=np.diff(B[SpkT[nInd]-PreCut-1:SpkT[nInd]-PreCut+Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]/3.)
				z[:,2:]-=np.tile(bA,(z.shape[0],1))#remove global changes
				#fine tuning based on correlations with individual channels
				bA-=np.mean(bA)
				z[goodCh,2:]-=np.outer((Vsbias[goodChInd]),bA2)#remove global changes II
				z[goodCh,1]-=Vsbias[goodChInd]*bA2[PreCut]#should change that in the detection instead...no! local vs. global measure...
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
				Cz=np.cumsum(z[:,PreCut-Lspike/2-1:PreCut+(Lspike+9)/2],axis=1)
				Pos=np.argmin(np.sum(Cz[:,Lspike:]+Cz[:,Lspike-1:-1]-Cz[:,:-Lspike]-Cz[:,1:-Lspike+1],axis=0))
				SAmp=-(Cz[:,Pos+Lspike]+Cz[:,Pos+Lspike-1]-Cz[:,Pos]-Cz[:,Pos+1])#/(LSpike-1)#pos. values
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
				if (iiii==0):#reverse clipping of central channel
					Oamp[0]=SAmp[0]#want that noisy channels are spatially limited (but would also bias results...)
				###remove 50 percentile (i.e. 37.5 percentile from former baseline)
				ShAmp0=np.sum(SAmp0)
				ShAmp=np.sum(Oamp)
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
	f.create_dataset('PeakIndex', data=PreCut)
	f.close()
	return

### Spike sorting
def readShapesFileInserted(TxtFile, HdfFile, NSpk, Sampling, Lspike=4):
	'''
	Iweights = np.array([[1000, 209, 272, 272, 272, 272, 209, 209, 209,0],\
	[750, 220, 271, 300, 271, 250, 199, 220, 199,0],\
	[600, 230, 266, 333, 266, 230, 189, 230, 189,0],\
	[500, 241, 259, 375, 259, 214, 180, 241, 180,0],\
	[750, 220, 300, 271, 250, 271, 220, 199, 199,0],\
	[679, 232, 297, 297, 248, 248, 208, 208, 190,0],\
	[572, 245, 291, 330, 245, 229, 197, 217, 182,0],\
	[486, 258, 282, 370, 240, 213, 187, 225, 174,0],\
	[600, 230, 333, 266, 230, 266, 230, 189, 189,0],\
	[572, 245, 330, 291, 229, 245, 217, 197, 182,0],\
	[514, 261, 321, 321, 227, 227, 204, 204, 175,0],\
	[454, 277, 309, 357, 223, 211, 193, 211, 167,0],\
	[500, 241, 375, 259, 214, 259, 241, 180, 180,0],\
	[486, 258, 370, 282, 213, 240, 225, 187, 174,0],\
	[454, 277, 357, 309, 211, 223, 211, 193, 167,0],\
	[414, 297, 339, 339, 208, 208, 198, 198, 161,0]])
	#Inn = np.array([0, 65, 64, 1, -64, -1, 63, -63, -65])
	Inn = np.argsort(np.array([4, 8, 7, 5, 1, 3, 6, 2, 0]))
	Ann = np.array([12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19])
	Iseq = np.array([65, 833, 1601, 2369, 3137, 3905, 577, 1345, 2145, 2913, 3681, 353, 1121, 1889, 2657, 3425, 69, 837, 1605, 2373, 3141, 3909, 581, 1349, 2149, 2917, 3685, 357, 1125, 1893, 2661, 3429, 73, 841, 1609, 2377, 3145, 3913, 585, 1353, 2153, 2921, 3689, 361, 1129, 1897, 2665, 3433, 77, 845, 1613, 2381, 3149, 3917, 589, 1357, 2157, 2925, 3693, 365, 1133, 1901, 2669, 3437, 81, 849, 1617, 2385, 3153, 3921, 593, 1361, 2161, 2929, 3697, 369, 1137, 1905, 2673, 3441, 85, 853, 1621, 2389, 3157, 3925, 597, 1365, 2165, 2933, 3701, 373, 1141, 1909, 2677, 3445, 89, 857, 1625, 2393, 3161, 3929, 601, 1369, 2169, 2937, 3705, 377, 1145, 1913, 2681, 3449, 93, 861, 1629, 2397, 3165, 3933, 605, 1373, 2173, 2941, 3709, 381, 1149, 1917, 2685, 3453, 97, 865, 1633, 2401, 3169, 3937, 609, 1377, 2113, 2881, 3649, 321, 1089, 1857, 2625, 3393, 101, 869, 1637, 2405, 3173, 3941, 613, 1381, 2117, 2885, 3653, 325, 1093, 1861, 2629, 3397, 105, 873, 1641, 2409, 3177, 3945, 617, 1385, 2121, 2889, 3657, 329, 1097, 1865, 2633, 3401, 109, 877, 1645, 2413, 3181, 3949, 621, 1389, 2125, 2893, 3661, 333, 1101, 1869, 2637, 3405, 113, 881, 1649, 2417, 3185, 3953, 625, 1393, 2129, 2897, 3665, 337, 1105, 1873, 2641, 3409, 117, 885, 1653, 2421, 3189, 3957, 629, 1397, 2133, 2901, 3669, 341, 1109, 1877, 2645, 3413, 121, 889, 1657, 2425, 3193, 3961, 633, 1401, 2137, 2905, 3673, 345, 1113, 1881, 2649, 3417, 125, 893, 1661, 2429, 3197, 3965, 637, 1405, 2141, 2909, 3677, 349, 1117, 1885, 2653, 3421])
	Template = np.array([[-406, 25101, 32871, 18504, 3709, -4793, -8225, -8875, -8292, -7287, -6233, -5275, -4456, -3777, -3220, -2767, -2397, -2095, -1848, -1643, -1474, -1332, -1214, -1113],\
	[-146, 26608, 32335, 17446, 2987, -5130, -8324, -8861, -8234, -7219, -6168, -5218, -4409, -3738, -3188, -2741, -2376, -2078, -1833, -1631, -1464, -1324, -1207, -1108],\
	[529, 28003, 31743, 16406, 2293, -5450, -8420, -8849, -8180, -7155, -6108, -5166, -4365, -3702, -3159, -2717, -2357, -2062, -1820, -1621, -1455, -1317, -1200, -1102],\
	[1554, 29279, 31091, 15386, 1626, -5756, -8511, -8838, -8130, -7096, -6051, -5117, -4324, -3668, -3132, -2695, -2339, -2048, -1809, -1611, -1447, -1310, -1195, -1098],\
	[2867, 30427, 30382, 14382, 985, -6047, -8595, -8826, -8082, -7039, -5998, -5070, -4285, -3636, -3106, -2674, -2322, -2034, -1798, -1602, -1440, -1304, -1190, -1094],\
	[4416, 31434, 29616, 13392, 369, -6320, -8672, -8812, -8033, -6983, -5945, -5024, -4247, -3605, -3081, -2654, -2306, -2021, -1787, -1594, -1433, -1298, -1185, -1090],\
	[6150, 32291, 28792, 12415, -223, -6578, -8740, -8793, -7982, -6925, -5891, -4978, -4209, -3574, -3056, -2634, -2290, -2008, -1776, -1585, -1426, -1293, -1181, -1086],\
	[8023, 32994, 27911, 11452, -791, -6818, -8798, -8768, -7928, -6866, -5837, -4932, -4170, -3543, -3030, -2613, -2273, -1994, -1765, -1576, -1418, -1286, -1176, -1082],\
	[9990, 33538, 26978, 10503, -1334, -7039, -8844, -8736, -7869, -6803, -5780, -4883, -4130, -3510, -3004, -2591, -2255, -1980, -1753, -1566, -1410, -1280, -1170, -1077],\
	[12009, 33924, 25997, 9570, -1852, -7243, -8878, -8697, -7805, -6737, -5720, -4832, -4088, -3475, -2976, -2569, -2236, -1965, -1741, -1556, -1401, -1272, -1164, -1072],\
	[14042, 34158, 24973, 8656, -2345, -7427, -8901, -8650, -7737, -6668, -5658, -4779, -4044, -3439, -2946, -2545, -2217, -1948, -1727, -1544, -1392, -1265, -1157, -1066],\
	[16060, 34241, 23921, 7764, -2811, -7594, -8913, -8598, -7664, -6596, -5594, -4725, -3999, -3403, -2916, -2520, -2196, -1932, -1713, -1533, -1382, -1256, -1150, -1060],\
	[18032, 34187, 22846, 6896, -3254, -7745, -8916, -8540, -7589, -6522, -5528, -4670, -3953, -3365, -2885, -2495, -2176, -1914, -1699, -1521, -1372, -1247, -1142, -1053],\
	[19936, 34009, 21757, 6055, -3671, -7881, -8912, -8478, -7512, -6448, -5463, -4614, -3907, -3327, -2854, -2469, -2155, -1897, -1684, -1508, -1362, -1239, -1135, -1047],\
	[21756, 33720, 20665, 5243, -4066, -8004, -8901, -8415, -7435, -6373, -5397, -4559, -3862, -3290, -2824, -2444, -2134, -1880, -1670, -1496, -1351, -1230, -1127, -1040],\
	[23479, 33335, 19576, 4460, -4439, -8118, -8888, -8351, -7359, -6301, -5334, -4506, -3818, -3254, -2794, -2420, -2114, -1863, -1656, -1484, -1341, -1221, -1120, -1034]])
	
	Iweights = np.array([[1000, 150, 200, 200, 200, 200, 150, 150, 150,0],\
	[666, 158, 198, 222, 198, 181, 142, 158, 142,0],\
	[500, 166, 195, 250, 195, 166, 135, 166, 135,0],\
	[400, 174, 189, 285, 189, 153, 128, 174, 128,0],\
	[666, 158, 222, 198, 181, 198, 158, 142, 142,0],\
	[585, 168, 220, 220, 180, 180, 149, 149, 135,0],\
	[472, 178, 215, 247, 178, 165, 140, 156, 129,0],\
	[387, 188, 207, 281, 174, 153, 132, 162, 123,0],\
	[500, 166, 250, 195, 166, 195, 166, 135, 135,0],\
	[472, 178, 247, 215, 165, 178, 156, 140, 129,0],\
	[414, 190, 240, 240, 163, 163, 146, 146, 123,0],\
	[356, 203, 229, 270, 160, 151, 137, 151, 118,0],\
	[400, 174, 285, 189, 153, 189, 174, 128, 128,0],\
	[387, 188, 281, 207, 153, 174, 162, 132, 123,0],\
	[356, 203, 270, 229, 151, 160, 151, 137, 118,0],\
	[320, 220, 255, 255, 149, 149, 142, 142, 113,0]])
	#Inn = np.array([0, 65, 64, 1, -64, -1, 63, -63, -65])
	Inn = np.argsort(np.array([4, 8, 7, 5, 1, 3, 6, 2, 0]))
	Ann = np.array([12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19])
	Iseq = np.array([65, 833, 1601, 2369, 3137, 3905, 577, 1345, 2145, 2913, 3681, 353, 1121, 1889, 2657, 3425, 69, 837, 1605, 2373, 3141, 3909, 581, 1349, 2149, 2917, 3685, 357, 1125, 1893, 2661, 3429, 73, 841, 1609, 2377, 3145, 3913, 585, 1353, 2153, 2921, 3689, 361, 1129, 1897, 2665, 3433, 77, 845, 1613, 2381, 3149, 3917, 589, 1357, 2157, 2925, 3693, 365, 1133, 1901, 2669, 3437, 81, 849, 1617, 2385, 3153, 3921, 593, 1361, 2161, 2929, 3697, 369, 1137, 1905, 2673, 3441, 85, 853, 1621, 2389, 3157, 3925, 597, 1365, 2165, 2933, 3701, 373, 1141, 1909, 2677, 3445, 89, 857, 1625, 2393, 3161, 3929, 601, 1369, 2169, 2937, 3705, 377, 1145, 1913, 2681, 3449, 93, 861, 1629, 2397, 3165, 3933, 605, 1373, 2173, 2941, 3709, 381, 1149, 1917, 2685, 3453, 97, 865, 1633, 2401, 3169, 3937, 609, 1377, 2113, 2881, 3649, 321, 1089, 1857, 2625, 3393, 101, 869, 1637, 2405, 3173, 3941, 613, 1381, 2117, 2885, 3653, 325, 1093, 1861, 2629, 3397, 105, 873, 1641, 2409, 3177, 3945, 617, 1385, 2121, 2889, 3657, 329, 1097, 1865, 2633, 3401, 109, 877, 1645, 2413, 3181, 3949, 621, 1389, 2125, 2893, 3661, 333, 1101, 1869, 2637, 3405, 113, 881, 1649, 2417, 3185, 3953, 625, 1393, 2129, 2897, 3665, 337, 1105, 1873, 2641, 3409, 117, 885, 1653, 2421, 3189, 3957, 629, 1397, 2133, 2901, 3669, 341, 1109, 1877, 2645, 3413, 121, 889, 1657, 2425, 3193, 3961, 633, 1401, 2137, 2905, 3673, 345, 1113, 1881, 2649, 3417, 125, 893, 1661, 2429, 3197, 3965, 637, 1405, 2141, 2909, 3677, 349, 1117, 1885, 2653, 3421])
	Template = np.array([[-406, 25101, 32871, 18504, 3709, -4793, -8225, -8875, -8292, -7287, -6233, -5275, -4456, -3777, -3220, -2767, -2397, -2095, -1848, -1643, -1474, -1332, -1214, -1113],\
	[-146, 26608, 32335, 17446, 2987, -5130, -8324, -8861, -8234, -7219, -6168, -5218, -4409, -3738, -3188, -2741, -2376, -2078, -1833, -1631, -1464, -1324, -1207, -1108],\
	[529, 28003, 31743, 16406, 2293, -5450, -8420, -8849, -8180, -7155, -6108, -5166, -4365, -3702, -3159, -2717, -2357, -2062, -1820, -1621, -1455, -1317, -1200, -1102],\
	[1554, 29279, 31091, 15386, 1626, -5756, -8511, -8838, -8130, -7096, -6051, -5117, -4324, -3668, -3132, -2695, -2339, -2048, -1809, -1611, -1447, -1310, -1195, -1098],\
	[2867, 30427, 30382, 14382, 985, -6047, -8595, -8826, -8082, -7039, -5998, -5070, -4285, -3636, -3106, -2674, -2322, -2034, -1798, -1602, -1440, -1304, -1190, -1094],\
	[4416, 31434, 29616, 13392, 369, -6320, -8672, -8812, -8033, -6983, -5945, -5024, -4247, -3605, -3081, -2654, -2306, -2021, -1787, -1594, -1433, -1298, -1185, -1090],\
	[6150, 32291, 28792, 12415, -223, -6578, -8740, -8793, -7982, -6925, -5891, -4978, -4209, -3574, -3056, -2634, -2290, -2008, -1776, -1585, -1426, -1293, -1181, -1086],\
	[8023, 32994, 27911, 11452, -791, -6818, -8798, -8768, -7928, -6866, -5837, -4932, -4170, -3543, -3030, -2613, -2273, -1994, -1765, -1576, -1418, -1286, -1176, -1082],\
	[9990, 33538, 26978, 10503, -1334, -7039, -8844, -8736, -7869, -6803, -5780, -4883, -4130, -3510, -3004, -2591, -2255, -1980, -1753, -1566, -1410, -1280, -1170, -1077],\
	[12009, 33924, 25997, 9570, -1852, -7243, -8878, -8697, -7805, -6737, -5720, -4832, -4088, -3475, -2976, -2569, -2236, -1965, -1741, -1556, -1401, -1272, -1164, -1072],\
	[14042, 34158, 24973, 8656, -2345, -7427, -8901, -8650, -7737, -6668, -5658, -4779, -4044, -3439, -2946, -2545, -2217, -1948, -1727, -1544, -1392, -1265, -1157, -1066],\
	[16060, 34241, 23921, 7764, -2811, -7594, -8913, -8598, -7664, -6596, -5594, -4725, -3999, -3403, -2916, -2520, -2196, -1932, -1713, -1533, -1382, -1256, -1150, -1060],\
	[18032, 34187, 22846, 6896, -3254, -7745, -8916, -8540, -7589, -6522, -5528, -4670, -3953, -3365, -2885, -2495, -2176, -1914, -1699, -1521, -1372, -1247, -1142, -1053],\
	[19936, 34009, 21757, 6055, -3671, -7881, -8912, -8478, -7512, -6448, -5463, -4614, -3907, -3327, -2854, -2469, -2155, -1897, -1684, -1508, -1362, -1239, -1135, -1047],\
	[21756, 33720, 20665, 5243, -4066, -8004, -8901, -8415, -7435, -6373, -5397, -4559, -3862, -3290, -2824, -2444, -2134, -1880, -1670, -1496, -1351, -1230, -1127, -1040],\
	[23479, 33335, 19576, 4460, -4439, -8118, -8888, -8351, -7359, -6301, -5334, -4506, -3818, -3254, -2794, -2420, -2114, -1863, -1656, -1484, -1341, -1221, -1120, -1034]])
	'''
	#h
	Iweights = np.array([[1000, 190, 250, 250, 250, 250, 190, 190, 190,0],\
	[1000, 203, 250, 285, 250, 222, 178, 203, 178,0],\
	[1000, 217, 250, 333, 250, 200, 166, 217, 166,0],\
	[666, 229, 247, 400, 247, 181, 156, 229, 156,0],\
	[1000, 203, 285, 250, 222, 250, 203, 178, 178,0],\
	[1000, 220, 285, 285, 222, 222, 188, 188, 168,0],\
	[1000, 238, 285, 333, 222, 200, 174, 198, 158,0],\
	[666, 255, 281, 400, 220, 181, 162, 207, 149,0],\
	[1000, 217, 333, 250, 200, 250, 217, 166, 166,0],\
	[1000, 238, 333, 285, 200, 222, 198, 174, 158,0],\
	[1000, 261, 333, 333, 200, 200, 182, 182, 150,0],\
	[666, 285, 326, 400, 198, 181, 168, 189, 142,0],\
	[666, 229, 400, 247, 181, 247, 229, 156, 156,0],\
	[666, 255, 400, 281, 181, 220, 207, 162, 149,0],\
	[666, 285, 400, 326, 181, 198, 189, 168, 142,0],\
	[585, 320, 387, 387, 180, 180, 174, 174, 135,0]])
	Inn = np.argsort(np.array([4, 8, 7, 5, 1, 3, 6, 2, 0]))
	#int[] Inn = new int[9] {0, 65, 64, 1, -64, -1, 63, -63, -65};
	#Ann = np.array([18, 17, 19, 16, 20, 15, 21, 14, 22, 13, 23, 12, 24, 11, 25])
	Ann = np.array([12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19])
	Iseq = np.array([65, 833, 1601, 2369, 3137, 3905, 577, 1345, 2145, 2913, 3681, 353, 1121, 1889, 2657, 3425, 69, 837, 1605, 2373, 3141, 3909, 581, 1349, 2149, 2917, 3685, 357, 1125, 1893, 2661, 3429, 73, 841, 1609, 2377, 3145, 3913, 585, 1353, 2153, 2921, 3689, 361, 1129, 1897, 2665, 3433, 77, 845, 1613, 2381, 3149, 3917, 589, 1357, 2157, 2925, 3693, 365, 1133, 1901, 2669, 3437, 81, 849, 1617, 2385, 3153, 3921, 593, 1361, 2161, 2929, 3697, 369, 1137, 1905, 2673, 3441, 85, 853, 1621, 2389, 3157, 3925, 597, 1365, 2165, 2933, 3701, 373, 1141, 1909, 2677, 3445, 89, 857, 1625, 2393, 3161, 3929, 601, 1369, 2169, 2937, 3705, 377, 1145, 1913, 2681, 3449, 93, 861, 1629, 2397, 3165, 3933, 605, 1373, 2173, 2941, 3709, 381, 1149, 1917, 2685, 3453, 97, 865, 1633, 2401, 3169, 3937, 609, 1377, 2113, 2881, 3649, 321, 1089, 1857, 2625, 3393, 101, 869, 1637, 2405, 3173, 3941, 613, 1381, 2117, 2885, 3653, 325, 1093, 1861, 2629, 3397, 105, 873, 1641, 2409, 3177, 3945, 617, 1385, 2121, 2889, 3657, 329, 1097, 1865, 2633, 3401, 109, 877, 1645, 2413, 3181, 3949, 621, 1389, 2125, 2893, 3661, 333, 1101, 1869, 2637, 3405, 113, 881, 1649, 2417, 3185, 3953, 625, 1393, 2129, 2897, 3665, 337, 1105, 1873, 2641, 3409, 117, 885, 1653, 2421, 3189, 3957, 629, 1397, 2133, 2901, 3669, 341, 1109, 1877, 2645, 3413, 121, 889, 1657, 2425, 3193, 3961, 633, 1401, 2137, 2905, 3673, 345, 1113, 1881, 2649, 3417, 125, 893, 1661, 2429, 3197, 3965, 637, 1405, 2141, 2909, 3677, 349, 1117, 1885, 2653, 3421])
	Template = np.array([[-406, 25101, 32871, 18504, 3709, -4793, -8225, -8875, -8292, -7287, -6233, -5275, -4456, -3777, -3220, -2767, -2397, -2095, -1848, -1643, -1474, -1332, -1214, -1113],\
	[-146, 26608, 32335, 17446, 2987, -5130, -8324, -8861, -8234, -7219, -6168, -5218, -4409, -3738, -3188, -2741, -2376, -2078, -1833, -1631, -1464, -1324, -1207, -1108],\
	[529, 28003, 31743, 16406, 2293, -5450, -8420, -8849, -8180, -7155, -6108, -5166, -4365, -3702, -3159, -2717, -2357, -2062, -1820, -1621, -1455, -1317, -1200, -1102],\
	[1554, 29279, 31091, 15386, 1626, -5756, -8511, -8838, -8130, -7096, -6051, -5117, -4324, -3668, -3132, -2695, -2339, -2048, -1809, -1611, -1447, -1310, -1195, -1098],\
	[2867, 30427, 30382, 14382, 985, -6047, -8595, -8826, -8082, -7039, -5998, -5070, -4285, -3636, -3106, -2674, -2322, -2034, -1798, -1602, -1440, -1304, -1190, -1094],\
	[4416, 31434, 29616, 13392, 369, -6320, -8672, -8812, -8033, -6983, -5945, -5024, -4247, -3605, -3081, -2654, -2306, -2021, -1787, -1594, -1433, -1298, -1185, -1090],\
	[6150, 32291, 28792, 12415, -223, -6578, -8740, -8793, -7982, -6925, -5891, -4978, -4209, -3574, -3056, -2634, -2290, -2008, -1776, -1585, -1426, -1293, -1181, -1086],\
	[8023, 32994, 27911, 11452, -791, -6818, -8798, -8768, -7928, -6866, -5837, -4932, -4170, -3543, -3030, -2613, -2273, -1994, -1765, -1576, -1418, -1286, -1176, -1082],\
	[9990, 33538, 26978, 10503, -1334, -7039, -8844, -8736, -7869, -6803, -5780, -4883, -4130, -3510, -3004, -2591, -2255, -1980, -1753, -1566, -1410, -1280, -1170, -1077],\
	[12009, 33924, 25997, 9570, -1852, -7243, -8878, -8697, -7805, -6737, -5720, -4832, -4088, -3475, -2976, -2569, -2236, -1965, -1741, -1556, -1401, -1272, -1164, -1072],\
	[14042, 34158, 24973, 8656, -2345, -7427, -8901, -8650, -7737, -6668, -5658, -4779, -4044, -3439, -2946, -2545, -2217, -1948, -1727, -1544, -1392, -1265, -1157, -1066],\
	[16060, 34241, 23921, 7764, -2811, -7594, -8913, -8598, -7664, -6596, -5594, -4725, -3999, -3403, -2916, -2520, -2196, -1932, -1713, -1533, -1382, -1256, -1150, -1060],\
	[18032, 34187, 22846, 6896, -3254, -7745, -8916, -8540, -7589, -6522, -5528, -4670, -3953, -3365, -2885, -2495, -2176, -1914, -1699, -1521, -1372, -1247, -1142, -1053],\
	[19936, 34009, 21757, 6055, -3671, -7881, -8912, -8478, -7512, -6448, -5463, -4614, -3907, -3327, -2854, -2469, -2155, -1897, -1684, -1508, -1362, -1239, -1135, -1047],\
	[21756, 33720, 20665, 5243, -4066, -8004, -8901, -8415, -7435, -6373, -5397, -4559, -3862, -3290, -2824, -2444, -2134, -1880, -1670, -1496, -1351, -1230, -1127, -1040],\
	[23479, 33335, 19576, 4460, -4439, -8118, -8888, -8351, -7359, -6301, -5334, -4506, -3818, -3254, -2794, -2420, -2114, -1863, -1656, -1484, -1341, -1221, -1120, -1034]])
	
	Template256=np.zeros((16,256),dtype=int)
	for i in range(16):
		Template256[i,np.arange(24)]=Template[i,:]
	
	#x[Iseq]=range
	#RIseq=np.argsort(Iseq)#(((Iseq/64)/16)*16+(Iseq/4)%16))#0...256, ranks, position to timelag%256, gives offset, need to subtract
	Rxseq=(np.arange(4096)/64/4)*16+((np.arange(4096)/4)%16)
	RIseq=np.argsort(Rxseq[Iseq])
	#print Rxseq
	#Ncut=int(Sampling)/1002+int(Sampling)/835 +7 #length of the cutouts
	#NcutL=int(Sampling)/501+int(Sampling)/835 +7 #length of the long cutouts
	#PreCut=6#5#define where the peak should be (needed for subtraction of median voltage)
	#Lspike=5#max. length of the spike
	n=0
	f=h5py.File(HdfFile,'r+')
	Sampling=int(f['Sampling'].value)
	Ncut=int(f['NCut'].value)
	NcutL=int(f['NCutLong'].value)
	PreCut=int(f['PreCut'].value)
	PostCut=int(f['PostCut'].value)
	Reverse=int(f['reverse Detection'].value)
	if Reverse:
		Bsupport=(2*PreCut)/3+PostCut/3+Ncut/6
		BInd=np.zeros((2,Bsupport),dtype=int)#maybe a bit long...
		BInd[:,:Ncut/6]=1#online baseline estimate
		BInd[:,Ncut/6:Ncut/6+(2*PreCut)/3]=2+np.arange((2*PreCut)/3,dtype=int)[None,:]
		BInd[0,-PostCut/3:]=NcutL+2+np.arange(-PostCut/3,0,dtype=int)
		BInd[1,-PostCut/3:]=Ncut+2+np.arange(-PostCut/3,0,dtype=int)
	else:
		Bsupport=PreCut/3+(2*PostCut)/3+Ncut/6
		BInd=np.zeros((2,Bsupport),dtype=int)#maybe a bit long...
		BInd[:,:Ncut/6]=1#online baseline estimate
		BInd[:,Ncut/6:Ncut/6+PreCut/3]=2+np.arange(PreCut/3,dtype=int)[None,:]
		BInd[0,-(2*PostCut)/3:]=NcutL+2+np.arange(-(2*PostCut)/3,0,dtype=int)
		BInd[1,-(2*PostCut)/3:]=Ncut+2+np.arange(-(2*PostCut)/3,0,dtype=int)
	#Lmax=2*int(Sampling)/1002#interval to look for sudden voltage jumps
	#list of shape histograms
	Slist=[]
	for i in range(240):
		Slist.append(np.zeros((512,Ncut-3,9),dtype=int))
	g=f['RawEvents']
	g.create_dataset('Locations', (NSpk,2))
	g.create_dataset('ShAmp0', (NSpk,))
	g.create_dataset('ShAmp', (NSpk,))
	g.create_dataset('ShAmpX', (NSpk,))
	g.create_dataset('ShArea', (NSpk,))
	g.create_dataset('goodCh', (NSpk,),dtype=bool)
	g.create_dataset('badCh', (NSpk,),dtype=int)
	g.create_dataset('Shapes', (NSpk,NcutL), fillvalue=0)
	#IgnSpk=np.array(g['PreSelectedEvents'].value,dtype=bool)
	IgnCh=np.array(f['NoisyChannels'].value,dtype=int)
	SortInd=np.array(g['SortInd'].value,dtype=int)
	SpkT=np.array(g['Times'].value)
	Qdavg=np.array(f['ChannelVariability'].value)
	recCh=np.array(f['recordedChannels'].value)
	
	PeakL=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	h=f['GlobalVoltageFluctuations']
	B=np.concatenate((3*h['medianVoltage'].value,np.zeros(100,dtype=int)+6141))
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
	#need a 9x16 matrix (9ch)
	#and 12x16 (4ch)
	I9=np.zeros((16,16,9),dtype=int)
	I12=np.zeros((16,16,12),dtype=int)
	H=np.zeros((7,7),dtype=int)+9
	H[1:-3,1:-3]=np.reshape(Inn,(3,3))
	#print H
	for ij1 in range(4):
		for ij2 in range(4):
				I9[ij1*4+ij2,:,:]=Iweights[:,H[A9[1,:]+1+ij1,A9[0,:]+1+ij2]]
				I12[ij1*4+ij2,:,:]=Iweights[:,H[A12[1,:]+1+ij1,A12[0,:]+1+ij2]]
	#print I9, I12
	fName=np.array(['_Shapes.txt','_ShapesX.txt'])
	'''
	BInd=np.zeros((2,int(Sampling)/1002+int(Sampling)/835+2),dtype=int)#maybe a bit long...
	BInd[:,:int(Sampling)/1002-3]=1
	BInd[:,int(Sampling)/1002-3:int(Sampling)/1002]=np.arange(2,5,dtype=int)[None,:]
	BInd[0,-(int(Sampling)/835)-2:]=np.arange(NcutL-(int(Sampling)/835),NcutL+2,dtype=int)
	BInd[1,-(int(Sampling)/835)-2:]=np.arange(Ncut-(int(Sampling)/835),Ncut+2,dtype=int)
	'''
	for iiii in range(2):
		b=file(TxtFile + fName[iiii])
		for i in b:
			z=np.array(i.split(),dtype=float)
			if not ((iiii==0)*(int(z[0]) in IgnCh)):
				Iind=((int(z[0])/64)%4)*4+(int(z[0])%4)
				#print Iind
				ICh=Rxseq[int(z[0])]#((int(z[0])/64)/16)*16+(int(z[0])/4)%16#need for Amplitude, i.e. correct time
				#print ICh, RIseq[ICh]
				#Iind=RIseq[ICh]%16
				nInd=SortIndInv[n]
				IT=int(SpkT[nInd])-3
				#print IT, IT%256, Iseq[IT%256],z[0], RIseq[ICh], ICh
				Itime=(+RIseq[ICh]-IT)%256#temporal offset (from where to start)
				Atime=(+(Itime+128)%256-128+IT)%15
				Itrange=np.arange(Itime-2,Itime+Ncut-PreCut+1)%256
				z=1.*np.reshape(z,(-1,Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])+2))
				#bad stuff #ignore channels in list
				goodCh=(np.min(z[:,2:]%4048,axis=1)>=48)#working channels
				badCh=np.nonzero((np.min(z[:,2:]%4048,axis=1)<48)+1*(np.in1d(z[:,0],IgnCh))+1*(z[:,0]==-1))[0]#not working or nonexistent channels
				goodChInd=np.array(z[goodCh,0],dtype=int)#indices for determining global fluctuations
				noCh=np.nonzero(z[:,0]==-1)[0]#channels that do not exist
				if iiii:
					z[:12,PreCut-1:Ncut+2]-=((Template256[RIseq[ICh]%16,Itrange])[None,:]\
					*(I12[Iind,RIseq[ICh]/16,:])[:,None])/Ann[Atime]/1000./64.
				else:
					z[:9,PreCut-1:Ncut+2]-=((Template256[RIseq[ICh]%16,Itrange])[None,:]\
					*(I9[Iind,RIseq[ICh]/16,:])[:,None])/Ann[Atime]/1000./64.
					#print ((Template256[RIseq[ICh]%16,Itrange])[None,:]\
					#*(I9[Iind,RIseq[ICh]/16,:])[:,None])/Ann[Atime]/1000./64.
				z[:,1]/=64.# was scaled version of 33 percentile of voltage
				###global fluctuations to subtract
				bA=B[SpkT[nInd]-PreCut:SpkT[nInd]-PreCut+Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]/3.
				bA2=np.diff(B[SpkT[nInd]-PreCut-1:SpkT[nInd]-PreCut+Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]/3.)
				z[:,2:]-=np.tile(bA,(z.shape[0],1))#remove global changes
				#fine tuning based on correlations with individual channels
				#bA-=np.mean(bA)
				z[goodCh,2:]-=np.outer((Vsbias[goodChInd]),bA2)#remove global changes II
				z[goodCh,1]-=Vsbias[goodChInd]*bA2[PreCut]#should change that in the detection instead...no! local vs. global measure...
				if (Iind==5) and not iiii:
					for ij in range(9):
						if goodCh[ij]:
							Slist[RIseq[ICh]/16*15+Atime][np.clip(np.array(z[ij,5:Ncut+2]*4,dtype=int)+400,0,511)\
							,np.arange(Ncut-3),ij]+=1
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
				Cz=np.cumsum(z[:,PreCut-Lspike/2-1:PreCut+(Lspike+9)/2],axis=1)
				Pos=np.argmin(np.sum(Cz[:,Lspike:]+Cz[:,Lspike-1:-1]-Cz[:,:-Lspike]-Cz[:,1:-Lspike+1],axis=0))
				SAmp=-(Cz[:,Pos+Lspike]+Cz[:,Pos+Lspike-1]-Cz[:,Pos]-Cz[:,Pos+1])#/(LSpike-1)#pos. values
				#look for sudden jumps
				#should not do that at all... would induce bias in location, which is worse than bias in amplitude
				#SAmp-=np.max(np.abs(np.diff(z[:,2:PreCut-Lspike/2+Pos],axis=1)),axis=1)\
				#+np.max(np.abs(np.diff(z[:,PreCut+(Lspike+1)/2+Pos+1:Lmax+1],axis=1)),axis=1)
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
				if (iiii==0):#reverse clipping of central channel
					Oamp[0]=SAmp[0]#want that noisy channels are spatially limited (but would also bias results...)
				###remove 50 percentile (i.e. 37.5 percentile from former baseline)
				ShAmp0=np.sum(SAmp0)
				ShAmp=np.sum(Oamp)
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
				if np.sum(1.*(np.array(z[goodCh,0]%4,dtype=int)==1)*(((np.array(z[goodCh,0],dtype=int)/64)%4)==1)):
					g['goodCh'][nInd]=True
				else:
					g['goodCh'][nInd]=False
				g['badCh'][nInd]=5-1*iiii-np.sum(1*goodCh[:5-1*iiii])
				g['Shapes'][nInd,:Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]=wA
				n+=1
		b.close()
	MedianShape=np.zeros((15,16,Ncut-3,9))
	MedianShapeDev=np.zeros((15,16,Ncut-3,9))
	Q=((np.arange(Ncut-3,dtype=int)-2)%256)
	TemplateMean=np.mean(Template256[:,Q],axis=0)
	for i in range(240):
		a=np.sum(Slist[i][:,0,:],axis=0)
		for j in range(9):
			MedianShape[i%15,i/15,:,j]=(np.argmax(np.diff(np.concatenate((np.zeros((1,Ncut-3))\
			,1*(np.cumsum(Slist[i][:,:,j],axis=0)>a[j]/2)),axis=0),axis=0),axis=0)-400)/4.
			MedianShapeDev[i%15,i/15,:,j]=-MedianShape[i%15,i/15,:,j]-(TemplateMean\
			*I9[5,i/15,j])*1./Ann[i%15]/1000./64.
	g['Locations'][:,:]=np.clip(g['Locations'].value\
	+1e-2*scipy.randn(g['Locations'].value.shape[0],2),0.,63.999)#to avoid discretization effects
	f.create_dataset('lenCutouts', data=Ncut)
	f.create_dataset('lenCutoutsLong', data=NcutL)
	f.create_dataset('PeakIndex', data=PreCut)
	f.create_dataset('EstimatedShapes', data=MedianShape)
	f.create_dataset('EstimatedShapeDev', data=MedianShapeDev)
	f.close()
	return

def MergeForwardBackward(HdfFileForward, HdfFileBackward, HdfOutput, IncludeLongSpikesForward=False, IncludeLongSpikesBackward=True, DFrames=3, MaxDist=1.):
	#should do both isolation of events and check whether there is a smaller event of the reverse detection
	#in the recpective spatio-temporal neighborhood.
	#shall output a new .hdf file (of smaller size)
	#maybe isolation before?
	MaxDist*=MaxDist
	f=h5py.File(HdfFileForward,'r+')
	g=f['RawEvents']
	Sampling=f['Sampling'].value
	PreCut=f['PreCut'].value
	PostCut=f['PostCut'].value
	NCut=f['NCut'].value
	tMax=f['tMax'].value
	nFrames=f['nFrames'].value
	ChannelVariability1=f['ChannelVariability'].value
	LocF=g['Locations'].value
	ShAmpX=g['ShAmpX'].value
	ShArea=g['ShArea'].value
	TimesF=g['Times'].value
	RepolarizingSpikesF=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	ShapesF=np.array(g['Shapes'].value)
	AmplitudesF=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/12.
	f.close()
	#
	f=h5py.File(HdfFileBackward,'r+')
	g=f['RawEvents']
	assert Sampling==f['Sampling'].value
	ChannelVariability2=f['ChannelVariability'].value
	LocB=g['Locations'].value
	ShAmpX=g['ShAmpX'].value
	ShArea=g['ShArea'].value
	TimesB=g['Times'].value
	RepolarizingSpikesB=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	ShapesB=np.array(g['Shapes'].value)
	AmplitudesB=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/12.
	f.close()
	#need single arrays of Times (sorted) and Amplitudes
	Times=np.concatenate((TimesF,TimesB))
	Ind=np.argsort(Times)
	LenF=TimesF.shape[0]
	FromF=(Ind<LenF)
	Amplitudes=np.concatenate((AmplitudesF,AmplitudesB))[Ind]
	RepolarizingSpikes=np.concatenate((RepolarizingSpikesF,RepolarizingSpikesB))[Ind]
	RepolarizingSpikes0=RepolarizingSpikes.copy()
	if IncludeLongSpikesForward:
		RepolarizingSpikes0=np.concatenate((np.ones(LenF,dtype=bool),RepolarizingSpikesB))[Ind]
	if IncludeLongSpikesBackward:
		RepolarizingSpikes0=np.concatenate((RepolarizingSpikesF,np.ones(TimesB.shape[0],dtype=bool)))[Ind]
	s=min(ShapesF.shape[1], ShapesB.shape[1])
	Shapes=np.concatenate((ShapesF[:,:s],ShapesB[:,:s]),axis=0)[Ind,:]
	Loc=np.concatenate((LocF,LocB),axis=0)[Ind,:]
	Times=Times[Ind]
	X=np.zeros((1001,2))
	Xa=np.zeros((1001))
	Xt=np.zeros((1001))
	Xf=np.zeros((1001),dtype=bool)
	if IncludeLongSpikesForward and IncludeLongSpikesBackward:
		Ispikes=np.ones(len(Times),dtype=bool)#isolated
		Dspikes=np.zeros(len(Times),dtype=bool)#both forward and reverse
		for i in range(500):
			X[i,:]=Loc[i,:]
			Xa[i]=Amplitudes[i]
			Xt[i]=int(Times[i])
			Xf[i]=FromF[i]
		for i in range(500,len(Times)):
			X[i%1001,:]=Loc[i,:]
			Xa[i%1001]=Amplitudes[i]
			Xt[i%1001]=int(Times[i])
			Xf[i%1001]=FromF[i]
			j=(i-500)%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
			Ind=np.nonzero((np.sum((X-X[j,:][None,:])**2,axis=1)<=MaxDist)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any((Xa[Ind]-Xa[j])>0):
				Ispikes[i-500]=False
			elif any(((Xa[Ind]-Xa[j])<=0)*(Xf[Ind]<>Xf[j])):
				Dspikes[i-500]=True
		for i in range(len(Times)-500,len(Times)):
			j=i%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
			Ind=np.nonzero((np.sum((X-X[j,:][None,:])**2,axis=1)<=MaxDist)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any((Xa[Ind]-Xa[j])>0):
				Ispikes[i]=False
			elif any(((Xa[Ind]-Xa[j])<=0)*(Xf[Ind]<>Xf[j])):
				Dspikes[i]=True
	else:
		Dspikes=np.zeros(len(Times),dtype=bool)
		Ispikes=RepolarizingSpikes.copy()
		RepolarizingInd=np.nonzero(RepolarizingSpikes0)[0]
		for i in range(500):
			k=RepolarizingInd[i]
			X[i,:]=Loc[k,:]
			Xa[i]=Amplitudes[k]
			Xt[i]=int(Times[k])
			Xf[i]=FromF[k]
		for i in range(500,np.sum(RepolarizingSpikes0)):
			k=RepolarizingInd[i]
			X[i%1001,:]=Loc[k,:]
			Xa[i%1001]=Amplitudes[k]
			Xt[i%1001]=int(Times[k])
			Xf[i%1001]=FromF[k]
			j=(i-500)%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
			Ind=np.nonzero((np.sum((X-X[j,:][None,:])**2,axis=1)<=MaxDist)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any((Xa[Ind]-Xa[j])>0):
				Ispikes[RepolarizingInd[i-500]]=False
			elif any(((Xa[Ind]-Xa[j])<=0)*(Xf[Ind]<>Xf[j])):
				Dspikes[RepolarizingInd[i-500]]=True
		for i in range(np.sum(RepolarizingSpikes0)-500,np.sum(RepolarizingSpikes0)):
			k=RepolarizingInd[i]
			j=i%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
			Ind=np.nonzero((np.sum((X-X[j,:][None,:])**2,axis=1)<=MaxDist)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any((Xa[Ind]-Xa[j])>0):
				Ispikes[k]=False
			elif any(((Xa[Ind]-Xa[j])<=0)*(Xf[Ind]<>Xf[j])):
				Dspikes[k]=True
	f=h5py.File(HdfOutput,'w')
	f.create_dataset('Sampling', data=Sampling)
	f.create_dataset('PreCut', data=PreCut)
	f.create_dataset('PostCut', data=PostCut)
	f.create_dataset('NCut', data=NCut)
	f.create_dataset('NCutLong', data=NCut)
	f.create_dataset('tMax', data=tMax)
	f.create_dataset('nFrames', data=nFrames)
	f.create_dataset('ChannelVariability1', data=ChannelVariability1)
	f.create_dataset('ChannelVariability2', data=ChannelVariability2)
	f.create_dataset('ChannelVariability', data=(ChannelVariability1+ChannelVariability2)/2.)
	f.create_dataset('Locations', data=Loc[Ispikes,:])
	f.create_dataset('Amplitudes', data=Amplitudes[Ispikes])
	f.create_dataset('Times', data=Times[Ispikes])
	f.create_dataset('RepolarizingSpikes', data=RepolarizingSpikes[Ispikes])
	f.create_dataset('Shapes', data=Shapes[Ispikes])
	f.create_dataset('Origin', data=FromF[Ispikes])
	f.create_dataset('Double', data=Dspikes[Ispikes])
	if IncludeLongSpikesForward or IncludeLongSpikesBackward:
		f.create_dataset('IncludeLongSpikes', data=True, dtype=bool)
	else:
		f.create_dataset('IncludeLongSpikes', data=False, dtype=bool)
	f.close()
	return
	

def IsolatedSpikes(HdfFile, IncludeLongSpikes=True, DFrames=3, MaxDist=1.):
	MaxDist*=MaxDist
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
	Amplitudes=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/12.#+1e-6*scipy.rand(Times.shape())#not like that!
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
			j=(i-500)%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
			Ind=np.nonzero(((Xa-Xa[j])>0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind,:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Ispikes[i-500]=False
		for i in range(len(Times)-500,len(Times)):
			j=i%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
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
			j=(i-500)%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
			Ind=np.nonzero(((Xa-Xa[j])>0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind,:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Ispikes[RepolarizingInd[i-500]]=False
		for i in range(np.sum(RepolarizingSpikes)-500,np.sum(RepolarizingSpikes)):
			k=RepolarizingInd[i]
			j=i%1001
			Ind=np.nonzero(((Xa-Xa[j])==0)*(np.abs(Xt-Xt[j])<=DFrames))[0]
			if any(np.sum((X[Ind[Ind<>j],:]-X[j,:][None,:])**2,axis=1)<=MaxDist):
				Xa[j]+=1e-6-2e-6*(scipy.rand(1)<0.5)
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
	g=f['RawEvents']
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
				NS=np.clip(np.sum(SpikesSI),20,1e12)
				#if NS<30:
				#	print NS
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
				pPoisson=np.clip(pPoisson*fNoise[i,ii]*NS/np.clip(np.sum(pPoisson*SpikesSiX),1e-12,1e12),0,1)#does not go well with clipping...
				pPoisson=np.clip(pPoisson*fNoise[i,ii]*NS/np.clip(np.sum(pPoisson*SpikesSiX),1e-12,1e12),0,1)#do a round of adjustment to fNoise
				pPoisson=np.clip(pPoisson*fNoise[i,ii]*NS/np.clip(np.sum(pPoisson*SpikesSiX),1e-12,1e12),0,1)#do a round of adjustment to fNoise
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

def MakeOldFormat(HdfFile,OldHdfFile):
	g=h5py.File(HdfFile,'r')
	f=h5py.File(OldHdfFile,'w')
	h=g['GlobalVoltageFluctuations']
	h0=g['RawEvents']
	h1=g['CorrelationAnalysis']
	h2=g['Cluster']
	f.create_dataset('Sampling', data=g['Sampling'].value)
	f.create_dataset('tMax', data=g['tMax'].value)
	f.create_dataset('nFrames', data=g['nFrames'].value)
	f.create_dataset('AmplitudeThresholds', data=g['AmplitudeThresholds'].value)
	f.create_dataset('RepolarizationThreshold', data=g['RepolarizationThreshold'].value)
	f.create_dataset('ignoreRecalibrations', data=g['ignoreRecalibrations'].value)
	f.create_dataset('lenCutouts', data=g['lenCutouts'].value)
	#f.create_dataset('NSpikes', data=)#doesn't make sense
	f.create_dataset('medianVoltage', data=h['medianVoltage'].value)
	f.create_dataset('recordedChannels', data=g['recordedChannels'].value)
	f.create_dataset('RecalibrationEvents', data=g['RecalibrationEvents'].value)
	f.create_dataset('Amplitudes', data=g['Amplitudes'].value)
	f.create_dataset('Times', data=g['Times'].value)
	f.create_dataset('Shapes', data=g['Shapes'].value)
	f.create_dataset('Locations', data=g['Locations'].value)
	f.create_dataset('Probability', data=h1['Probability'].value)
	Ispikes=h0['IsolatedSpikes'].value
	f.create_dataset('Channels', data=(h0['Channels'].value)[Ispikes])
	f.create_dataset('ShAmp', data=(h0['ShAmpX'].value)[Ispikes])
	f.create_dataset('ShArea', data=(h0['ShArea'].value)[Ispikes])
	
	f0=f.create_group('CorrelationAnalysis')
	f0.create_dataset('Amplitudes', data=g['Amplitudes'].value)
	f0.create_dataset('Cluster', data=h1['Cluster'].value)
	#Map=np.repeat(np.repeat(np.reshape(np.arange(4096),(64,64)),2,axis=1),2,axis=0).flatten()
	#f0.create_dataset('Channels', data=Map[h1['Units'].value])#different dimension-->map
	f0.create_dataset('Channels', data=h2['ClusterId'].value)#different dimension-->from clustering...
	f0.create_dataset('Times', data=g['Times'].value)
	f0.create_dataset('Corr', data=g['Corr'].value)
	f0.create_dataset('Probability', data=h1['Probability'].value)
	f0.create_dataset('Noise', data=h2['Noise'].value.flatten())#different dimension
	#f0.create_dataset('pValue', data=g['pValue'].value.flatten())#different dimension-->doesn't make sense
	
	f1=f.create_group('Cluster')
	f1.create_dataset('ClusterId', data=h2['ClusterId'].value)
	f1.create_dataset('NCount', data=h2['NCount'].value)
	f1.create_dataset('ClusterVariance', data=h2['ClusterVariance'].value)
	f1.create_dataset('CLocations', data=h2['CLocations'].value)
	f1.create_dataset('gradientMap', data=h2['gradientMap'].value)
	f1.create_dataset('CAreaMatrix', data=h2['CAreaMatrix'].value)
	f1.create_dataset('CDensityDecay', data=h2['CDensityDecay'].value)
	f1.create_dataset('CBoundaries', data=h2['CBoundaries'].value)
	h2h=h2['Parameter']
	fh=f.create_group('Parameter')
	fh.create_dataset('nBins', data=h2h['nBins'])
	fh.create_dataset('nNext', data=h2h['fNext']*g['tMax'].value)
	fh.create_dataset('nNextMin', data=h2h['fNextMin'].value*g['tMax'].value)
	fh.create_dataset('gradientThreshold', data=h2h['gradientThreshold'].value)
	g.close()
	f.close()
	return

def MakeSmallHdf(HdfFile,SmallHdfFile):
	g=h5py.File(HdfFile,'r')
	f=h5py.File(SmallHdfFile,'w')
	f.create_dataset('Sampling', data=g['Sampling'].value)
	f.create_dataset('tMax', data=g['tMax'].value)
	f.create_dataset('nFrames', data=g['nFrames'].value)
	f.create_dataset('AmplitudeThresholds', data=g['AmplitudeThresholds'].value)
	f.create_dataset('RepolarizationThreshold', data=g['RepolarizationThreshold'].value)
	f.create_dataset('ignoreRecalibrations', data=g['ignoreRecalibrations'].value)
	f.create_dataset('lenCutouts', data=g['lenCutouts'].value)
	f.create_dataset('recordedChannels', data=g['recordedChannels'].value, compression='lzf')
	f.create_dataset('Amplitudes', data=g['Amplitudes'].value, compression='lzf')
	f.create_dataset('Times', data=g['Times'].value, compression='lzf')
	f.create_dataset('Locations', data=g['Locations'].value, compression='lzf')
	f.create_dataset('RepolarizingSpikes', data=g['RepolarizingSpikes'].value, compression='lzf')
	f.create_dataset('IncludeLongSpikes', data=g['IncludeLongSpikes'].value)
	f.create_dataset('RecalibrationOffsets', data=g['RecalibrationOffsets'].value, compression='lzf')
	g.close()
	f.close()
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
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=(g['Amplitudes'].value)[SingleSpk]
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

def TempBias(HdfFile,NoisyChFile,Ns0=5./12, minAmp=2,Res=12,nNext=500, tau0=50., tauL=200.):
	f=h5py.File(NoisyChFile,'r+')
	NoisyAreas=f['NoisyAreas'].value
	f.close()
	nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=(g['Amplitudes'].value)[SingleSpk]
	Times=(g['Times'].value)[SingleSpk]
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res0=int(nBins0/64)
	nBins=64*Res
	Ns=int(Ns0*Res)
	tau=1000./Sampling*1./tau0
	taul=42./tauL
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
	Times=Times[(Mask[Spikes]<>-1)*(Amp>minAmp)]
	Spikes=Mask[Spikes[(Mask[Spikes]<>-1)*(Amp>minAmp)]]
	SpikesC=Loc[:,1]+1j*Loc[:,0]
	N=Spikes.shape[0]
	#make smoothened histograms of surrounding spikes
	PPB=[]
	for i in range(nNext,N-nNext):
		D=np.exp(-np.abs(SpikesC[i-nNext:i+nNext+1]-SpikesC[i])*taul)
		T=np.exp(-np.abs(Times[i-nNext:i+nNext+1]-Times[i])*tau)
		PPB.append((np.sum(D[:nNext]*T[:nNext])-np.sum(D[-nNext:]*T[-nNext:]))*1./(np.sum(D*T)-1))
	TB=np.histogram2d(Loc[nNext:N-nNext,0]*Res,Loc[nNext:N-nNext,1]*Res\
	,bins=(np.arange(nBins+1),np.arange(nBins+1)),weights=(np.array(PPB)+1.)/2.)[0]
	Nloc=np.histogram2d(Loc[nNext:-nNext,0]*nBins/64.,Loc[nNext:-nNext,1]*nBins/64.,\
	bins=(np.arange(nBins+1),np.arange(nBins+1)))[0]
	#smoothening
	H0=np.concatenate((np.zeros((1,nBins)),np.cumsum(TB,axis=0)))
	H0=H0[Ns:,:]-H0[:-Ns,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:]-H1[:,:-Ns]
	TB1=H1*1./Ns**2
	H0=np.concatenate((np.zeros((1,nBins)),np.cumsum(Nloc,axis=0)))
	H0=H0[Ns:,:]-H0[:-Ns,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:]-H1[:,:-Ns]
	Nloc1=H1*1./Ns**2
	if 'TempBias' in g:
		del g['TempBias']
	g.create_dataset('TempBias', data=(TB1+0.5*(Nloc1<1.))*1./(Nloc1+(Nloc1<1.)))#want to avoid noise
	g.close()
	return

#is surround firing rate higher at some locations?
def FRBias(HdfFile,NoisyChFile,Ns0=5./12, minAmp=2,Res=12,nNext=500, tau0=50., tauL=500.):
	f=h5py.File(NoisyChFile,'r+')
	NoisyAreas=f['NoisyAreas'].value
	f.close()
	nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=(g['Amplitudes'].value)[SingleSpk]
	Times=(g['Times'].value)[SingleSpk]
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res0=int(nBins0/64)
	nBins=64*Res
	Ns=int(Ns0*Res)
	tau=1000./Sampling*1./tau0
	taul=42./tauL
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
	Times=Times[(Mask[Spikes]<>-1)*(Amp>minAmp)]
	Spikes=Mask[Spikes[(Mask[Spikes]<>-1)*(Amp>minAmp)]]
	SpikesC=Loc[:,1]+1j*Loc[:,0]
	N=Spikes.shape[0]
	#make smoothened histograms of surrounding spikes
	PPB=[]
	for i in range(nNext,N-nNext):
		D=np.exp(-np.abs(SpikesC[i-nNext:i+nNext+1]-SpikesC[i])*taul)
		T=np.exp(-np.abs(Times[i-nNext:i+nNext+1]-Times[i])*tau)
		PPB.append(np.clip((np.sum(D*T)-1)*2./nNext,0,1))
	TB=np.histogram2d(Loc[nNext:N-nNext,0]*Res,Loc[nNext:N-nNext,1]*Res\
	,bins=(np.arange(nBins+1),np.arange(nBins+1)),weights=(np.array(PPB)))[0]
	Nloc=np.histogram2d(Loc[nNext:-nNext,0]*nBins/64.,Loc[nNext:-nNext,1]*nBins/64.,\
	bins=(np.arange(nBins+1),np.arange(nBins+1)))[0]
	#smoothening
	H0=np.concatenate((np.zeros((1,nBins)),np.cumsum(TB,axis=0)))
	H0=H0[Ns:,:]-H0[:-Ns,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:]-H1[:,:-Ns]
	TB1=H1*1./Ns**2
	H0=np.concatenate((np.zeros((1,nBins)),np.cumsum(Nloc,axis=0)))
	H0=H0[Ns:,:]-H0[:-Ns,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1)),np.cumsum(H0,axis=1)),axis=1)
	H1=H1[:,Ns:]-H1[:,:-Ns]
	Nloc1=H1*1./Ns**2
	if 'FRBias' in g:
		del g['FRBias']
	g.create_dataset('FRBias', data=(TB1+0.5*(Nloc1<1.))*1./(Nloc1+(Nloc1<1.)))#want to avoid noise
	g.close()
	return

def PeriEventActivity2(HdfFile,NoisyChFile, minAmp=2,Res=8,Ns0=3./8,Na0=4, nNext=200, Ares=50):
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
	Nsx=int(Ns0*Res)/2
	Nsy=int(Ns0*Res-1)/2
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
	#Dhist=np.zeros((N,5))
	#shall increase resolution here
	#directly make a histogram
	Aloc=np.zeros((nBins,nBins,Ares))
	Dloc=np.zeros((nBins,nBins,Ares))
	Nloc=np.zeros((nBins,nBins))
	for i in range(nNext,N-nNext):
		y=int(Loc[i,0]*nBins/64.)
		x=int(Loc[i,1]*nBins/64.)
		D=np.abs(SpikesC[i-nNext:i+nNext+1]-SpikesC[i])
		Dx=(D<27)*(D>2)
		Dl1=np.cumsum(Dx)
		Dl1-=Dl1[nNext]
		Ind=np.nonzero(np.abs(Dl1<31)*Dx)[0]
		Dc=np.clip(((D[Ind]-2.)*2.),0,Ares-1)
		Nc=len(Dc)
		Dloc[y,x,:]+=np.histogram((Dc)%Ares,bins=np.arange(Ares+1),weights=1.-Dc%1)[0]*1./(Nc+(Nc==0))
		Dloc[y,x,:]+=np.histogram((Dc+1)%Ares,bins=np.arange(Ares+1),weights=Dc%1)[0]*1./(Nc+(Nc==0))
		Ac=np.angle(SpikesC[i-nNext:i+nNext+1][Ind]-SpikesC[i])*Ares*0.5/np.pi-0.5
		Aloc[y,x,:]+=np.histogram(Ac%Ares,bins=np.arange(Ares+1),weights=1.-Ac%1)[0]*1./(Nc+(Nc==0))
		Aloc[y,x,:]+=np.histogram((Ac+1)%Ares,bins=np.arange(Ares+1),weights=Ac%1)[0]*1./(Nc+(Nc==0))
		Nloc[y,x]+=1
	#smoothening
	if Nsy and Nsx:
		X=np.concatenate((np.arange(Nsx)[::-1],np.arange(nBins),np.arange(nBins-Nsy+1,nBins)[::-1]))
	elif Nsx:
		X=np.concatenate((np.arange(Nsx)[::-1],np.arange(nBins)))
	else:
		X=np.arange(nBins)
	H0=np.concatenate((np.zeros((1,nBins,5)),np.cumsum(Dloc[X,:,:],axis=0)))
	H0=H0[Ns:,:,:]-H0[:-Ns,:,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1,5)),np.cumsum(H0[:,X,:],axis=1)),axis=1)
	H1=H1[:,Ns:,:]-H1[:,:-Ns,:]
	Dloc1=H1*1./Ns**2
	H0=np.concatenate((np.zeros((1,nBins,6)),np.cumsum(Aloc[X,:,:],axis=0)))
	H0=H0[Ns:,:,:]-H0[:-Ns,:,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1,6)),np.cumsum(H0[:,X,:],axis=1)),axis=1)
	H1=H1[:,Ns:,:]-H1[:,:-Ns,:]
	Aloc1=H1*1./Ns**2
	H0=np.concatenate((np.zeros((1,nBins)),np.cumsum(Nloc[X,:],axis=0)))
	H0=H0[Ns:,:]-H0[:-Ns,:]
	H1=np.concatenate((np.zeros((nBins-Ns+1,1)),np.cumsum(H0[:,X],axis=1)),axis=1)
	H1=H1[:,Ns:]-H1[:,:-Ns]
	Nloc1=H1*1./Ns**2
	#Normalization: maybe not needed...
	Dloc1=Dloc1*1./np.clip(Nloc1,1e-6,1e12)[:,:,None]
	Aloc1=Aloc1*1./np.clip(Nloc1,1e-6,1e12)[:,:,None]
	#local averages
	X=np.concatenate((np.arange(Na+1)[::-1],np.arange(nBins),np.arange(nBins-Na,nBins)[::-1]))
	H0=np.cumsum((np.sum(Dloc1[X,:,:],axis=2)>0)*1.,axis=0)
	H0=H0[2*Na+1:,:]-H0[:-2*Na-1,:]
	H1=np.cumsum(H0[:,X],axis=1)
	Navg=np.clip(H1[:,2*Na+1:]-H1[:,:-2*Na-1],1,1e12)
	H0=np.cumsum(Dloc1[X,:,:],axis=0)
	H0=H0[2*Na+1:,:,:]-H0[:-2*Na-1,:,:]
	H1=np.cumsum(H0[:,X,:],axis=1)
	H1=H1[:,2*Na+1:,:]-H1[:,:-2*Na-1,:]
	Davg=H1*1./Navg[:,:,None]
	H0=np.cumsum(Aloc1[X,:,:],axis=0)
	H0=H0[2*Na+1:,:,:]-H0[:-2*Na-1,:,:]
	H1=np.cumsum(H0[:,X,:],axis=1)
	H1=H1[:,2*Na+1:,:]-H1[:,:-2*Na-1,:]
	Aavg=H1*1./Navg[:,:,None]
	#need to find a normalized average to compare with (instead of second recording)
	#renormalize averaged Histogram
	Anorm=Aavg*Nloc1[:,:,None]*1./Navg[:,:,None]
	Dnorm=Davg*Nloc1[:,:,None]*1./Navg[:,:,None]
	#do that comparison. save data in a separate file?!
	#first Amplitudes, then distances...
	D=np.cumsum(np.concatenate(((Aloc1-Anorm)[:,:,:4][:,:,::-1]\
	,(Aloc1-Anorm),(Aloc1-Anorm)[:,:,-3:][:,:,::-1]),axis=2),axis=2)
	dD=D[:,:,5:]-D[:,:,:-5]
	N0=dD.shape[0]
	N1=dD.shape[1]
	N2=dD.shape[2]
	#DN=np.zeros((N0,N1,N2))
	DN=dD.copy()
	RdistP=dD[:,:,1:-1]
	RdistN=dD[:,:,1:-1]
	#RdistP=np.zeros((N0,N1,50))
	#increases in #spikes
	for i in range(3):
		if i:
			X=np.cumsum(np.concatenate((np.zeros((1,dD.shape[1],dD.shape[2])),dD),axis=0),axis=0)
			Y=np.cumsum(np.concatenate((np.zeros((dD.shape[0]-i,1,dD.shape[2])),X[(i+1):,:,:]-X[:-(i+1),:,:]),axis=1),axis=1)
			DN[(i)/2:N0-(i+1)/2,(i)/2:N1-(i+1)/2,:]=(Y[:,(i+1):,:]-Y[:,:-(i+1),:])*1./(i+1)
		Rdist=DN[:,:,1:-1]
		for k in range(i+1):
			for j in range(i+1):
				Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:]=np.min(np.concatenate(\
				(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:][:,:,:,None]),axis=3),axis=3)
		RdistP=np.min(np.concatenate((RdistP[:,:,:,None],Rdist[:,:,:,None]),axis=3),axis=3)
		Rdist=DN[:,:,1:-1]
		for k in range(i+1):
			for j in range(i+1):
				Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:]=np.max(np.concatenate(\
				(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:][:,:,:,None]),axis=3),axis=3)
		RdistN=np.max(np.concatenate((RdistN[:,:,:,None],Rdist[:,:,:,None]),axis=3),axis=3)
	#need to make that positive
	RdistP=np.clip(RdistP,0,1e12)
	#need to make that negative
	RdistN=np.clip(RdistN,-1e12,0)
	#stats? --> don't have a resolved version of the number of spikes,
	A=np.cumsum(np.concatenate(((Aloc1)[:,:,:4][:,:,::-1],(Aloc1),(Aloc1)[:,:,-3:][:,:,::-1]),axis=2),axis=2)
	#N=np.zeros((N0,N1,N2))
	N=np.min(np.concatenate(((A[:,:,5:-2]-A[:,:,:-7])[:,:,:,None]\
	,(A[:,:,6:-1]-A[:,:,1:-6])[:,:,:,None],(A[:,:,7:]-A[:,:,2:-5])[:,:,:,None]),axis=3),axis=3)
	#save angles
	f=h5py.File(ExcFile,'r+')
	f.create_dataset('AngleP%i' %recId, data=RdistP)
	f.create_dataset('AngleN%i' %recId, data=RdistN)
	f.create_dataset('AngleNN%i' %recId, data=N)
	f.close()
	#distances
	D=np.cumsum(np.concatenate(((Dloc1-Dnorm)[:,:,:4][:,:,::-1]\
	,(Dloc1-Dnorm),(Dloc1-Dnorm)[:,:,-3:][:,:,::-1]),axis=2),axis=2)
	dD=D[:,:,5:]-D[:,:,:-5]
	N0=dD.shape[0]
	N1=dD.shape[1]
	N2=dD.shape[2]
	#DN=np.zeros((N0,N1,N2))
	DN=dD.copy()
	RdistP=dD[:,:,1:-1]
	RdistN=dD[:,:,1:-1]
	#RdistP=np.zeros((N0,N1,50))
	#increases in #spikes
	for i in range(3):
		if i:
			X=np.cumsum(np.concatenate((np.zeros((1,dD.shape[1],dD.shape[2])),dD),axis=0),axis=0)
			Y=np.cumsum(np.concatenate((np.zeros((dD.shape[0]-i,1,dD.shape[2])),X[(i+1):,:,:]-X[:-(i+1),:,:]),axis=1),axis=1)
			DN[(i)/2:N0-(i+1)/2,(i)/2:N1-(i+1)/2,:]=(Y[:,(i+1):,:]-Y[:,:-(i+1),:])*1./(i+1)
		Rdist=DN[:,:,1:-1]
		for k in range(i+1):
			for j in range(i+1):
				Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:]=np.min(np.concatenate(\
				(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:][:,:,:,None]),axis=3),axis=3)
		RdistP=np.min(np.concatenate((RdistP[:,:,:,None],Rdist[:,:,:,None]),axis=3),axis=3)
		Rdist=DN[:,:,1:-1]
		for k in range(i+1):
			for j in range(i+1):
				Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:]=np.max(np.concatenate(\
				(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:][:,:,:,None]),axis=3),axis=3)
		RdistN=np.max(np.concatenate((RdistN[:,:,:,None],Rdist[:,:,:,None]),axis=3),axis=3)
	#need to make that positive
	RdistP=np.clip(RdistP,0,1e12)
	#need to make that negative
	RdistN=np.clip(RdistN,-1e12,0)
	#stats? --> don't have a resolved version of the number of spikes,
	A=np.cumsum(np.concatenate(((Dloc1)[:,:,:4][:,:,::-1],(Dloc1),(Dloc1)[:,:,-3:][:,:,::-1]),axis=2),axis=2)
	#N=np.zeros((N0,N1,N2))
	N=np.min(np.concatenate(((A[:,:,5:-2]-A[:,:,:-7])[:,:,:,None]\
	,(A[:,:,6:-1]-A[:,:,1:-6])[:,:,:,None],(A[:,:,7:]-A[:,:,2:-5])[:,:,:,None]),axis=3),axis=3)
	#save distances
	f=h5py.File(ExcFile,'r+')
	f.create_dataset('DistancesP%i' %recId, data=RdistP)
	f.create_dataset('DistancesN%i' %recId, data=RdistN)
	f.create_dataset('DistancesNN%i' %recId, data=N)
	f.close()
	return

#should not care about time and just look for the locations of surrounding spikes.
#will lead to some measure how spreadout the activity is --> silent periods > activity periods
#not very quantitative... and dependent on the location on the chip!!!

def Excitability(HdfFile,NoisyChFile, minAmp=2.5, tau0=0.02, NNext=10, minFreq=0.1):
	f=h5py.File(NoisyChFile,'r+')
	NoisyAreas=f['NoisyAreas'].value
	f.close()
	nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	f=g['SurrogateData']
	NewIndMap=f['NewIndMap'].value
	OldInd0=f['OldInd'].value
	Loc=g['Locations'].value[NewIndMap,:]
	Amp=g['Amplitudes'].value[NewIndMap]
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)[NewIndMap]
	Times=f['NewTimes'].value
	g.close()
	#print Loc.shape, Amp.shape, Times.shape
	Res0=int(nBins0/64)
	Spikes=np.clip(np.array(Res0*(Loc[:,1]),dtype=int),0,nBins0-1)\
	+np.clip(np.array(Res0*(Loc[:,0]),dtype=int),0,nBins0-1)*(nBins0)
	#Exclusion of small spikes and noisy channels
	Ind=(True-NoisyAreas[Spikes])*(Amp>minAmp)*SingleSpk
	OldInd1=np.zeros(Ind.shape,dtype=bool)
	OldInd1[OldInd0]=True#shall ignore surrogate spikes
	OldInd=np.nonzero(OldInd1[Ind])[0]
	#print OldInd0.shape, OldInd.shape
	Loc=Loc[Ind,:]
	#Amp=Amp[Ind]
	Spikes=Spikes[Ind]
	Times=Times[Ind]
	N=Spikes.shape[0]
	print N
	#Parameter
	tau=tau0*Sampling#frames#20ms
	minN=int(minFreq*tMax*4./Res0**2)
	#print minN
	
	NSpikes=np.histogram2d(Loc[OldInd,0]*Res0,Loc[OldInd,1]*Res0,bins=(np.arange(nBins0+1),np.arange(nBins0+1)))[0]
	#make smoothened histograms of surrounding spikes
	NSpikesS=np.zeros(NSpikes.shape, dtype=int)
	NSpikesS[1:,1:]=NSpikes[:-1,:-1]+NSpikes[:-1,1:]+NSpikes[1:,:-1]+NSpikes[1:,1:]
	#print NSpikes.max(), NSpikes.mean(), NSpikes.shape, NSpikes.sum()
	minN=max(np.percentile(NSpikesS,80),minN)
	RefInd=NSpikesS.flatten()>minN
	RefIndi=np.nonzero(RefInd)[0]
	RefIndLoc=np.argwhere(NSpikesS>minN)
	#print len(RefIndi)
	NRef=np.sum(RefInd)
	RefIndL=np.zeros(RefInd.shape,dtype=int)
	RefIndL[RefIndi]=np.arange(NRef)
	RefIndSrci=np.unique(np.array([RefIndi,RefIndi-nBins0-1,RefIndi-nBins0,RefIndi-1]))
	RefIndSrc=np.zeros(RefInd.shape,dtype=bool)
	RefIndSrc[RefIndSrci]=True
	NRefSrc=np.sum(RefIndSrc)
	RefIndSrcL=np.zeros(RefInd.shape,dtype=int)
	RefIndSrcL[RefIndSrci]=np.arange(NRefSrc)
	print time.localtime(time.time())[:6]
	#need a map
	RefMap=[]
	for i in range(NRefSrc):
		x=RefIndSrci[i]+np.array([0,nBins0+1,nBins0,1])
		RefMap.append(RefIndL[x[np.in1d(x,RefIndi)]])
	#need an activity list for each reference channel
	print len(RefMap)
	RefAList=[]
	for i in range(NRef):
		RefAList.append([])
	#print len(RefAList)
	#first round of percentile estimation
	OldIndX=OldInd[(OldInd>=500)*(OldInd<N-500)]
	for i in OldIndX[RefIndSrc[Spikes[OldIndX]]]:#OldInd
		x=np.array(Loc[i,:],dtype=int)
		y=np.array(np.delete(Loc[i-500:i+500,:],500,0))
		Lag=np.exp(-np.abs(Times[i]-np.delete(Times[i-500:i+500],500))*1./tau)#temporal
		q=np.sum(np.exp(-np.sqrt((y[:,0]-x[0])**2+(y[:,1]-x[1])**2)*2./NNext)*Lag)
		for k in RefMap[RefIndSrcL[Spikes[i]]]:
			RefAList[k].append(q)
	p=np.linspace(0,1,100)[1:]
	pMat=np.zeros((NRef,99))
	#compute percentiles
	for i in range(NRef):
		ql=len(RefAList[i])
		if ql:
			qs=np.sort(np.array(RefAList[i]))+1e-6*scipy.rand(ql)
			pMat[i,:]=np.interp(p,np.linspace(0,1.,ql),qs)
	print pMat.max(), pMat.shape
	#lists of surrounding reference channels
	RefList=[]
	DistList=[]
	DistH0=np.zeros((nBins0**2,NNext),dtype=int)
	for i in range(nBins0**2):
		Dist=np.sqrt(((RefIndi%nBins0)-0.5-(i%nBins0))**2+((RefIndi/nBins0)-0.5-(i/nBins0))**2)
		RefList.append(np.nonzero(Dist<2*NNext)[0])
		DistList.append(np.array(Dist[RefList[i]]/2,dtype=int))
		DistH0[i,:]=np.bincount(DistList[i],minlength=NNext)
	DistH=np.clip(DistH0,1,1e12)
	#for each spike...
	N=len(OldInd)
	SpikesE=np.zeros((N,NNext))
	#for i in OldIndX[RefIndSrc[Spikes[OldIndX]]]:
	for i in np.arange(N)[(OldInd>=500)*(OldInd<N-500)]:
		j=OldInd[i]
		x=np.array(Loc[j,:],dtype=int)
		y=np.array(np.delete(Loc[j-500:j+500,:],500,0))
		Lag=np.exp(-np.abs(Times[j]-np.delete(Times[j-500:j+500],500))*1./tau)#temporal
		q=np.sum(np.exp(-np.sqrt((y[:,0]-x[0])**2+(y[:,1]-x[1])**2)*2./NNext)*Lag)+1e-6*scipy.rand(1)
		pq=np.sum(pMat[RefList[Spikes[j]],:]<q,axis=1)+0.5
		SpikesE[i,:]=np.histogram(DistList[Spikes[j]],bins=np.arange(NNext+1),weights=pq)[0]*1./DistH[Spikes[j],:]
	for i in np.arange(N)[(OldInd>=N-500)]:#boundaries
		j=OldInd[i]
		x=np.array(Loc[j,:],dtype=int)
		y=np.array(np.delete(Loc[j-500:,:],500,0))
		Lag=np.exp(-np.abs(Times[j]-np.delete(Times[j-500:],500))*1./tau)#temporal
		q=np.sum(np.exp(-np.sqrt((y[:,0]-x[0])**2+(y[:,1]-x[1])**2)*2./NNext)*Lag)+1e-6*scipy.rand(1)
		pq=np.sum(pMat[RefList[Spikes[j]],:]<q,axis=1)+0.5
		SpikesE[i,:]=np.histogram(DistList[Spikes[j]],bins=np.arange(NNext+1),weights=pq)[0]*1./DistH[Spikes[j],:]
	for i in np.arange(N)[(OldInd<500)]:
		j=OldInd[i]
		x=np.array(Loc[j,:],dtype=int)
		y=np.array(np.delete(Loc[:j+500,:],j,0))
		Lag=np.exp(-np.abs(Times[j]-np.delete(Times[:j+500],j))*1./tau)#temporal
		q=np.sum(np.exp(-np.sqrt((y[:,0]-x[0])**2+(y[:,1]-x[1])**2)*2./NNext)*Lag)+1e-6*scipy.rand(1)
		pq=np.sum(pMat[RefList[Spikes[j]],:]<q,axis=1)+0.5
		SpikesE[i,:]=np.histogram(DistList[Spikes[j]],bins=np.arange(NNext+1),weights=pq)[0]*1./DistH[Spikes[j],:]
	#need to figure out distance dependent bias
	#Histogram for each location, normalize
	#average over locations (ignoring zeros ==> count, fraction, inverse multiplication)
	NS=np.histogram(Spikes[OldInd],bins=np.arange(nBins0**2+1))[0]
	Nzero=np.clip(np.sum((DistH0*(NS>0)[:,None])==0,axis=0),1,1e12)
	NS=np.clip(NS,1,1e12)
	Eavg=np.zeros(NNext)
	for i in range(NNext):
		Eavg[i]=np.sum((np.histogram(Spikes[OldInd],bins=np.arange(nBins0**2+1)\
		,weights=SpikesE[:,i])[0])*1./NS)*1./(nBins0**2-Nzero[i])
	print np.mean(Eavg)
	DistN=np.sum(DistH0<>0,axis=1)
	print np.mean(DistN)
	relE=np.sum(SpikesE*1./Eavg[None,:],axis=1)*1./np.clip(DistN[Spikes[OldInd]],1,1e12)
	print np.mean(SpikesE), np.mean(relE)
	relEAll=np.zeros(OldInd0.shape)#rather want something of standard shape...
	relEAll[Ind[OldInd0]]=relE
	g=h5py.File(HdfFile,'r+')
	if 'Excitability3' in g:
		del g['Excitability3']
	f=g.create_group('Excitability3')
	f.create_dataset('relE', data=relEAll)
	f.create_dataset('Indices', data=Ind[OldInd0])
	f.create_dataset('minAmp', data=minAmp)
	f.create_dataset('NNext', data=NNext)
	f.create_dataset('minFreq', data=minFreq)
	f.create_dataset('RefInd', data=RefInd)
	#f.create_dataset('Loc', data=Loc)
	#f.create_dataset('Times', data=Times)
	g.close()
	return

def SurrogateData(HdfFile,NoisyChFile, Nbs=10000, Nscale=10, nNext=50\
, minAmp=2, minFreq=0.1, aHistory=0.9):
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
	Times=g['Times'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	g.close()
	Res0=int(nBins0/64)
	Spikes=np.clip(np.array(Res0*(Loc[:,1]),dtype=int),0,nBins0-1)\
	+np.clip(np.array(Res0*(Loc[:,0]),dtype=int),0,nBins0-1)*(nBins0)
	#Exclusion of small spikes and noisy channels and long spikes
	Ind=(True-NoisyAreas[Spikes])*(Amp>minAmp)*SingleSpk
	#Loc=Loc[Ind,:]
	#Amp=Amp[Ind]
	#Spikes=Spikes[Ind]
	#Times=Times[Ind]
	nNext0=int(np.sqrt(nNext))
	N=np.sum(Ind)
	#Parameter
	minN=int(minFreq*tMax*4./Res0**2)
	#find reference units
	NSpikes=np.histogram2d(Loc[Ind,0]*Res0,Loc[Ind,1]*Res0,bins=(np.arange(nBins0+1),np.arange(nBins0+1)))[0]
	#make smoothened histograms of surrounding spikes
	NSpikesS=np.zeros(NSpikes.shape, dtype=int)
	NSpikesS[1:,1:]=NSpikes[:-1,:-1]+NSpikes[:-1,1:]+NSpikes[1:,:-1]+NSpikes[1:,1:]
	minN=max(np.percentile(NSpikesS,80),minN)
	RefInd=NSpikesS.flatten()>minN
	RefIndi=np.nonzero(RefInd)[0]
	#locations for clustering
	RefIndLoc=np.argwhere(NSpikesS>minN)
	#number of reference units
	NRef=np.sum(RefInd)
	print NRef
	print N
	print Ind.shape
	#clustering
	ClusterDist=scipy.cluster.hierarchy.distance.pdist(RefIndLoc,'euclidean')
	ClusterLinkage=scipy.cluster.hierarchy.linkage(ClusterDist,method='weighted')
	FCluster=scipy.cluster.hierarchy.fcluster(ClusterLinkage,500,'maxclust')
	#find periods of similar activity
	print time.localtime(time.time())[:6]
	SpikesL=np.zeros((2*Nbs,2*nNext),dtype=int)
	MapL=np.zeros((nBins0**2),dtype=int)-1
	MapL[RefInd]=FCluster
	Ind0=Ind.copy()
	#which spikes are reference spikes
	Ind[Ind]*=MapL[Spikes[Ind]]<>-1
	MaskS=np.nonzero(Ind)[0]
	#only reference spikes
	SpikesLmsk=MapL[Spikes[MaskS]]
	RelInd=np.zeros((len(SpikesLmsk),Nscale),dtype=int)
	RelInd[:,0]=np.arange(len(SpikesLmsk))
	#kernel for comparison
	Hweights=np.hamming(2*nNext0+1)#[nNext0+1:]
	#need a history dependence
	SpkOverlap=scipy.rand(2*Nbs,2)*0.001#do some jittering to avoid equal cases
	SpkHMask=np.zeros((2*Nbs))
	for i in range(2*nNext):
		SpikesL[:,i]=SpikesLmsk[i:i+2*Nbs]
	#print SpikesL
	print len(SpikesLmsk)
	#fixed SpikesL
	for i in range(nNext,Nbs):#(fix initial and final parts later)
		SpkOverlap[1:,i%2]=SpkOverlap[:-1,(i-1)%2]
		SpkOverlap[0,i%2]=0
		SpkOverlap*=aHistory#current data represent 10% of total...
		Spk0=np.sum((SpikesL\
		-SpikesLmsk[None,i-nNext:nNext+i])==0,axis=1)
		SpkOverlap[:,i%2]+=np.convolve(Spk0,Hweights,mode='same')
		'''
		for k in range(1,nNext0):
			SpkOverlap[:,i%2]+=Hweights[k]\
			*np.sum((SpikesL-SpikesLmsk[i+k-nNext:i+nNext+k][None,:])==0,axis=1)
			SpkOverlap[:,i%2]+=Hweights[k]\
			*np.sum((SpikesL-SpikesLmsk[i-nNext-k:i+nNext-k][None,:])==0,axis=1)
		'''
		
		#recursively remove smaller values to find local maxima
		SpkOverlapX=SpkOverlap[:,i%2]*(SpkHMask<=0)#remove previous maxima
		SpkOMask=np.ones(2*Nbs,dtype=bool)
		for k in np.arange(1,25):
			SpkOMask[k:]*=SpkOverlapX[:-k]<SpkOverlapX[k:]
			SpkOMask[:-k]*=SpkOverlapX[k:]<SpkOverlapX[:-k]
		SpkOMask*=SpkOverlap[:,(i+1)%2]<SpkOverlap[:,i%2]
		#boundaries
		SpkOMask[:2*nNext0]=False
		SpkOMask[-2*nNext0:]=False
		#identity
		SpkOMask[np.arange(i-nNext+1-nNext,i)%(2*Nbs)]=False
		#iInd=np.arange(i,i+2*nNext,dtype=int)%(2*Nbs)
		#iInd=(iInd+1)%(2*Nbs)
		#SpkOMask[iInd]=False
		RelInd[i,1:]=np.argsort(SpkOMask*SpkOverlap[:,i%2])[-Nscale+1:]+nNext
		#set points that have been used to zero... bias??? -- doesn't matter if there are more good maxima
		SpkHMask-=1
		SpkHMask[RelInd[i,1:]-nNext]=nNext#indices that have been used already
	print RelInd[i-20:i+1,:], i
	print time.localtime(time.time())[:6]
	#changing SpikesL
	for i in range(Nbs,len(SpikesLmsk)-Nbs-2*nNext):
		if ((i+Nbs)%(2*Nbs))>1:
			SpkOverlap[1:(i+Nbs)%(2*Nbs),i%2]=SpkOverlap[:(i+Nbs)%(2*Nbs)-1,(i-1)%2]
		if ((i+Nbs)%(2*Nbs))<(2*Nbs-1):
			SpkOverlap[(i+Nbs)%(2*Nbs)+1:,i%2]=SpkOverlap[(i+Nbs)%(2*Nbs):-1,(i-1)%2]
		SpkOverlap[(i+Nbs)%(2*Nbs),i%2]=0#used wrong history...
		#SpkOverlap[-1,i%2]=0#don't need as I'm updating SpikesL as well?--> do need as I'm not shifting!!!
		SpkOverlap*=aHistory#current data represent 10% of total...
		SpikesL[(i+Nbs)%(2*Nbs),:]=SpikesLmsk[i+Nbs:i+2*nNext+Nbs]
		#SpkOverlap[:,i%2]+=np.sum((SpikesL-SpikesLmsk[None,i-nNext:i+nNext])==0,axis=1)
		Spk0=np.sum((SpikesL-SpikesLmsk[None,i-nNext:nNext+i])==0,axis=1)
		SpkOverlap[:,i%2]+=np.convolve(Spk0,Hweights,mode='same')
		'''
		for k in range(1,nNext0):
			SpkOverlap[:,i%2]+=Hweights[k]*np.sum((SpikesL-SpikesLmsk[i+k-nNext:i+nNext+k][None,:])==0,axis=1)
			SpkOverlap[:,i%2]+=Hweights[k]*np.sum((SpikesL-SpikesLmsk[i-k-nNext:i+nNext-k][None,:])==0,axis=1)
		'''
		#recursively remove smaller values to find local maxima
		SpkOverlapX=SpkOverlap[:,i%2]*(SpkHMask<=0)
		SpkOMask=np.ones(2*Nbs,dtype=bool)
		for k in np.arange(1,25):
			SpkOMask[k:]*=SpkOverlapX[:-k]<SpkOverlapX[k:]
			SpkOMask[:-k]*=SpkOverlapX[k:]<SpkOverlapX[:-k]
			SpkOMask[:k]*=SpkOverlapX[-k:]<SpkOverlapX[:k]
			SpkOMask[-k:]*=SpkOverlapX[:k]<SpkOverlapX[-k:]
		SpkOMask*=SpkOverlap[:,(i+1)%2]<SpkOverlap[:,i%2]
		#boundaries
		SpkOMask[(np.arange(-2*nNext0,2*nNext0)+i+Nbs)%(2*Nbs)]=False
		#identity
		SpkOMask[np.arange(i-nNext+1-nNext,i)%(2*Nbs)]=False
		RelInd[i,1:]=np.argsort(SpkOMask*SpkOverlap[:,i%2])[-Nscale+1:]
		SpkHMask-=1
		SpkHMask[RelInd[i,1:]]=nNext
		RelInd[i,1:]=(RelInd[i,1:]-i+Nbs)%(2*Nbs)-Nbs+i+nNext
		#set points that have been used to zero... bias??? -- doesn't matter if there are more good maxima
	print RelInd[i-20:i+1,:], i
	print time.localtime(time.time())[:6]
	i0=len(SpikesLmsk)-Nbs-1-2*nNext
	#fixed SpikesL
	for i in range(i0+1,len(SpikesLmsk)-nNext):
		if ((i0+Nbs)%(2*Nbs))>1:
			SpkOverlap[1:(i0+Nbs)%(2*Nbs),i%2]=SpkOverlap[:(i0+Nbs)%(2*Nbs)-1,(i-1)%2]
		if ((i0+Nbs)%(2*Nbs))<(2*Nbs-1):
			SpkOverlap[(i0+Nbs)%(2*Nbs)+1:,i%2]=SpkOverlap[(i0+Nbs)%(2*Nbs):-1,(i-1)%2]
		SpkOverlap[(i0+Nbs)%(2*Nbs),i%2]=0
		SpkOverlap*=aHistory#current data represent 10% of total...
		#SpkOverlap[:,i%2]+=np.sum((SpikesL-SpikesLmsk[None,i-nNext:i+nNext])==0,axis=1)
		Spk0=np.sum((SpikesL-SpikesLmsk[None,i-nNext:nNext+i])==0,axis=1)
		SpkOverlap[:,i%2]+=np.convolve(Spk0,Hweights,mode='same')
		'''
		for k in range(1,nNext0):
			SpkOverlap[:,i%2]+=Hweights[k]*np.sum((SpikesL-SpikesLmsk[i+k-nNext:i+nNext+k][None,:])==0,axis=1)
			SpkOverlap[:,i%2]+=Hweights[k]*np.sum((SpikesL-SpikesLmsk[i-k-nNext:i+nNext-k][None,:])==0,axis=1)
		'''
		#recursively remove smaller values to find local maxima
		SpkOverlapX=SpkOverlap[:,i%2]*(SpkHMask<=0)
		SpkOMask=np.ones(2*Nbs,dtype=bool)
		for k in np.arange(1,25):
			SpkOMask[k:]*=SpkOverlapX[:-k]<SpkOverlapX[k:]
			SpkOMask[:-k]*=SpkOverlapX[k:]<SpkOverlapX[:-k]
			SpkOMask[:k]*=SpkOverlapX[-k:]<SpkOverlapX[:k]
			SpkOMask[-k:]*=SpkOverlapX[:k]<SpkOverlapX[-k:]
		SpkOMask*=SpkOverlap[:,(i+1)%2]<SpkOverlap[:,i%2]
		#boundaries
		SpkOMask[(np.arange(-2*nNext0,2*nNext0)+i0+Nbs)%(2*Nbs)]=False
		#identity
		SpkOMask[np.arange(i-nNext+1-nNext,i)%(2*Nbs)]=False
		SpkHMask-=1
		RelInd[i,1:]=np.argsort(SpkOMask*SpkOverlap[:,i%2])[-Nscale+1:]
		SpkHMask[RelInd[i,1:]]=nNext
		RelInd[i,1:]=(RelInd[i,1:]-i0+Nbs)%(2*Nbs)-Nbs+i0+nNext
		#set points that have been used to zero... bias??? -- doesn't matter if there are more good maxima
	print RelInd[i-20:i+1,:], i
	#final
	RelInd[:nNext,:]=(np.arange(nNext)[:,None]+RelInd[nNext,:][None,:])\
	-nNext
	RelInd[len(SpikesLmsk)-nNext:,:]=(np.arange(1,nNext+1)[:,None]\
	+RelInd[len(SpikesLmsk)-nNext-1,:][None,:])
	#construct a new (artificial) spike train
	#have all indices of reference spikes before a gap
	#make a list, jitter spikes from other places and sort
	#NewSpikes=[]
	#NewLoc0=[]
	#NewLoc1=[]
	NewTimes=[]
	#NewAmp=[]
	NewIndMap=[]
	for x in range(len(Times)):
		NewIndMap.append(x)
		#NewLoc0.append(Loc[x,0])
		#NewLoc1.append(Loc[x,1])
		#NewAmp.append(Amp[x])
		NewTimes.append(Times[x])
	for i in range(len(SpikesLmsk)):
		t0=Times[MaskS[RelInd[i,0]]]
		for k in range(1,Nscale):
			ti=Times[MaskS[RelInd[i,k]]]
			for x in range(MaskS[RelInd[i,k]],MaskS[RelInd[i,k]+1]):
				if Ind0[x]:#only use large spikes, not long and not from noisy channels.
					NewIndMap.append(x)
					#NewLoc0.append(Loc[x,0])
					#NewLoc1.append(Loc[x,1])
					#NewAmp.append(Amp[x])
					NewTimes.append(Times[x]+scipy.rand(1)-0.5+t0-ti)
	#N=len(NewSpikes)
	print len(NewTimes)
	NewIndMap=np.array(NewIndMap)
	NewTimes=np.array(NewTimes)
	#NewAmp=np.array(NewAmp)
	#NewLoc=np.concatenate((np.array(NewLoc0)[:,None],np.array(NewLoc1)[:,None]),axis=1)
	#sorting
	NewInd=np.argsort(NewTimes)
	OldInd=np.nonzero(NewInd<len(Times))[0]
	NewTimes=NewTimes[NewInd]
	NewIndMap=NewIndMap[NewInd]
	#NewLoc=NewLoc[NewInd,:]
	#NewAmp=NewAmp[NewInd]
	#Stats
	#Histogram of raw data usage
	H0=np.bincount(RelInd.flatten())
	H1=np.bincount(H0, minlength=1000)
	H1[999]=np.sum(H1[999:])
	H1=H1[:1000]
	#firing rate dependency
	ISI=Times[MaskS[2*nNext:]]-Times[MaskS[:-2*nNext]]
	ISI=np.concatenate((np.ones(nNext)*ISI[0],ISI,np.ones(nNext)*ISI[-1]))
	H2=np.histogram2d(np.clip(ISI,0,Sampling-1),np.clip(H0,0,999)\
	,bins=(np.arange(Sampling),np.arange(1000)))[0]
	g=h5py.File(HdfFile,'r+')
	if 'SurrogateData' in g:
		del g['SurrogateData']
	f=g.create_group('SurrogateData')
	f.create_dataset('nNext', data=nNext)
	f.create_dataset('minAmp', data=minAmp)
	f.create_dataset('minFreq', data=minFreq)
	f.create_dataset('Nscale', data=Nscale)
	f.create_dataset('aHistory', data=aHistory)
	f.create_dataset('FCluster', data=FCluster)
	f.create_dataset('RefInd', data=RefInd)
	#f.create_dataset('NewLoc', data=NewLoc)
	f.create_dataset('NewTimes', data=NewTimes)
	#f.create_dataset('NewAmplitudes', data=NewAmp)
	f.create_dataset('OldInd', data=OldInd)
	f.create_dataset('NewIndMap', data=NewIndMap)
	f.create_dataset('SourceHist', data=H1)
	f.create_dataset('RateDependency', data=H2)
	g.close()
	return
