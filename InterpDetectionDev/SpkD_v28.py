import numpy as np
import h5py
import scipy.stats
import scipy.cluster
import time
#from IPython import get_ipython
#ipython = get_ipython()

#ipython.run_line_magic(u'load_ext', u'cythonmagic')

###This is the functions used in the program that writes the .txt Files into a .hdf Files
# and estimates the spatial origins of spikes, marks events detected in multiple channels,
# runs a correlation analysis and, based on these results, does a clustering.

### read Info file
def readInfoFile(TxtFile, HdfFile):
	X=[]
	Y=[]
	recCh=[]
	#k=0
	#n=0
	n1=0
	b=file(TxtFile + '_Info.txt')
	VarName=''
	ThrScale=2
	Ascale=-64
	SqIv=0
	for i in b:
		if '#' in i:
			VarName=i
			#k+=1
			print i
			if VarName=='# Recording channels4:\n':
				XCh4=np.zeros((len(recCh),12),dtype=int)
				n0=0
			if VarName=='# Recording channels5:\n':
				XCh5=np.zeros((len(recCh),9),dtype=int)
				SIprod=np.zeros((len(recCh),13),dtype=long)
				QdAvg=np.zeros((len(recCh)),dtype=int)
				QdAvgN=np.zeros((len(recCh)),dtype=int)
				n0=0
		else:
			if VarName=='# Number of frames:\n':
				nFrames=int(i)
				print i
			if VarName=='# Duration (s):\n':
				tMax=float(i)
				print i
			if VarName=='# Sampling rate:\n':
				Sampling=float(i[:-1])
				print i
			if VarName=='# Threshold scaling:\n':
				ThrScale=int(i)
				print i
			if VarName=='# Amplitude scaling:\n':
				Ascale=int(i)
				print i
			if VarName=='# Smoothing window:\n':
				SmoothN=int(i)
				print i
			if VarName=='#Number of spikes (4 channel):\n':
				print i
			if VarName=='#Number of spikes (5 channel):\n':
				print i
			if VarName==('# Detection threshold*%i:\n'%ThrScale):
				SpkThreshold=np.array(i.split(),dtype=int)*1./ThrScale
				print i
			if VarName=='# Detection threshold:\n':
				SpkThreshold=np.array(i.split(),dtype=int)*1./ThrScale
				print i
			if VarName==('# Repolarization threshold*%i:\n'%ThrScale):
				SpkRPthreshold=int(i)*1./ThrScale
				print i
			if VarName=='# Repolarization threshold:\n':
				SpkRPthreshold=int(i)*1./ThrScale
				print i
			if VarName=='# Recalibration trigger:\n':
				recalibrationTrigger=int(i)
				print i
			if VarName=='# Cutouts:\n':#CutPre, CutPost, tCut, tCutLong, df  in frames
				cutoutDim=np.array(i.split(),dtype=int)
				print i
			if VarName=='# Recording channels:\n':
				recCh.append(int(i))
			if VarName=='# Recording channel4:\n':
				XCh4[n0,:]=np.array(i.split(),dtype=int)
				n0+=1
			if VarName=='# Recording channels5:\n':
				XCh5[n0,:]=np.array(i.split(),dtype=int)
				n0+=1
			if VarName=='# Recalibration events:\n':
				z=np.array(i.split(),dtype=int)
				X.append(z[0])
				QdAvg+=z[1:]
				QdAvgN+=z[1:]>0
				#n+=1
			if VarName=='#Sum(squared global fluctuations):\n':
				SqIglobal=np.array(i.split(),dtype=long)
			if VarName=='#Sum(squared channel fluctuations):\n':
				SqIv=np.array(i.split(),dtype=long)
			if VarName=='#Sum(product of channel and global fluctuations):	\n':
				SIprod[:,n1]=np.array(i.split(),dtype=long)
				n1+=1
			if VarName=='#Sum(avg. deviations from global fluctuations):\n':
				Vsbias=-np.array(i.split(),dtype=long)*1./Sampling/nFrames#in ADCcounts
			if VarName=='#Sum(avg. squared deviations from global fluctuations):\n':
				Vsqbias=np.array(i.split(),dtype=long)*1./Sampling**2/nFrames
	b.close()
	Recalib=np.array(X)# was used for recalibration events
	Qdavg=np.clip(QdAvg,10*np.abs(Ascale),1e12)*SmoothN*1./np.clip(QdAvgN,1,1e12)*1./np.abs(Ascale)#/n#in ADCcounts (for individual channels!)
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
	f.create_dataset('ThrScale', data=ThrScale)
	f.create_dataset('Ascale', data=Ascale/SmoothN)
	f.create_dataset('SmoothN', data=SmoothN)
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
def readAvgFile(TxtFile, HdfFile, NCh):
	b=file(TxtFile + '_Avg.txt')
	X=[]
	for i in b:
		X.append(int(i))
	b.close()
	B=np.array(X)
	Bstd=np.std(B)/4.
	f=h5py.File(HdfFile,'r+')
	nFrames=f['nFrames'].value
	#Ascale=int(f['Ascale'].value)
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
	if 'reverse Detection' in f:
		if f['reverse Detection'].value:
			g.create_dataset('medianVoltage', data=B[::-1]/2.)
		else:
			g.create_dataset('medianVoltage', data=B/2.)
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
	g=h5py.File(HdfFile,'r+')
	ThrScale=g['ThrScale'].value
	for i in b:
		z=np.array(i.split(),dtype=int)
		X.append(z[0])
		Y.append(z[1])
		Amp.append((z[2])*1./ThrScale)
		Y1.append(z[3])
		YR.append(z[4])
		#Y2.append(9+(np.sum((z[5:8]+z[6:])/2)+(z[5]+z[8])/2)*3)
		#Y2.append(z[5])
	b.close()
	b=file(TxtFile + '_SpikesX.txt')
	for i in b:
		z=np.array(i.split(),dtype=int)
		X.append(z[0]+NCh)
		Y.append(z[1])
		Amp.append((z[2])*1./ThrScale)
		Y1.append(z[3])
		YR.append(z[4])
		#Y2.append(12+(np.sum(z[5:9])))
		#Y2.append(z[5])
	b.close()
	SpkCh=np.array(X)
	SpkT=np.array(Y)+scipy.rand(SpkCh.shape[0])#want to randomize events that are in same frame.
	SpkAmp=np.array(Amp)
	PeakL=np.array(Y1,dtype=bool)
	#Ntraces=np.array(Y2)#not using that one here, maybe should store it anyway
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
	
	if 'RawEvents' in g:
		del g['RawEvents']
	f=g.create_group('RawEvents')
	f.create_dataset('SortInd', data=Ind, dtype=int)#to find spikes and raw data in the .txt files
	f.create_dataset('Amplitudes', data=SpkAmp[PreSelectedEvents][Ind])
	f.create_dataset('Channels', data=SpkCh[PreSelectedEvents][Ind])
	f.create_dataset('Times', data=SpkT[PreSelectedEvents][Ind])
	f.create_dataset('RepolarizingSpikes', data=PeakL[PreSelectedEvents][Ind] ,dtype=bool)
	f.create_dataset('RecalibrationOffsets', data=RecalibOffset[PreSelectedEvents][Ind])
	#f.create_dataset('NumberOfCutoutChannels', data=Ntraces[PreSelectedEvents][Ind])
	g.create_dataset('NoisyChannels', data=IgnCh, dtype=int)
	f.create_dataset('PreSelectedEvents', data=PreSelectedEvents, dtype=bool)
	g.close()
	return NSpk

#ipython.run_cell_magic(u'cython',u'',(u'cdef int n = 1234\n' 
#u'print(n)'))

class LocationFinder:
	def __init__(self,eCentered,Ascale,IgnCh,PreCut,Ncut,NcutL,Qdavg):
		self.Qdavg=Qdavg
		self.Vscale=2.
		self.PreCut=PreCut
		self.Ncut=Ncut
		self.Ncut0=Ncut+4
		self.NcutD=NcutL-self.Ncut
		self.NcutL0=NcutL+4
		self.Ascale=Ascale
		self.eCentered=eCentered
		self.UseCh=np.ones(4097,dtype=bool)
		self.UseCh[-1]=False
		self.UseCh[IgnCh]=False
		if eCentered:
			self.Mapx=[np.array([5,1,6,2,7,3,8,4]),np.array([6,2,0,4,5]),np.array([1,6,7,3,0]),\
			np.array([4,0,2,7,8]),np.array([5,1,0,3,8]),np.array([1,0,4]),\
			np.array([2,0,1]),np.array([0,2,3]),np.array([4,0,3])]
			self.MapxL=np.array([8,5,5,5,5,3,3,3,3])
			self.Ax=np.array([[0,0,1,0,-1,-1,1,1,-1],[0,-1,0,1,0,-1,-1,1,1]])
			self.Nch=9
			self.Nch0=5
			self.OampNorm=np.ones(self.Nch,dtype=float)
			self.OampNorm[1:]+=1.
			self.OampNorm[5:]+=1.
			self.Lx=5*self.Ncut0
			self.LxL=5*self.NcutL0
		else:
			self.Mapx=[np.array([4,5,1,2,3,10,11]),np.array([4,5,6,7,2,3,0]),np.array([0,1,3,6,7,8,9]),\
			np.array([0,1,2,8,9,10,11]),np.array([0,1,11,5]),np.array([0,1,4,6]),\
			np.array([1,2,5,7]),np.array([2,1,6,8]),np.array([2,3,7,9]),\
			np.array([2,3,8,10]),np.array([0,3,11,9]),np.array([10,0,3,4])]
			self.MapxL=np.array([7,7,7,7,4,4,4,4,4,4,4,4])
			self.Ax=np.array([[0,1,1,0,0,1,2,2,1,0,-1,-1],[0,0,1,1,-1,-1,0,1,2,2,1,0]])
			self.Nch=12
			self.Nch0=4
			self.OampNorm=2*np.ones(self.Nch,dtype=float)
			self.OampNorm[4:]+=1.
			self.Lx=4*self.Ncut0
			self.LxL=4*self.NcutL0
		self.Ib=np.arange(self.Nch,dtype=int)*4
		self.Ib[:self.Nch0]+=np.arange(self.Nch0,dtype=int)*Ncut
		self.Ib[self.Nch0:]+=self.Nch0*Ncut
		self.IbL=np.arange(self.Nch,dtype=int)*4
		self.IbL[:self.Nch0]+=np.arange(self.Nch0,dtype=int)*NcutL
		self.IbL[self.Nch0:]+=self.Nch0*NcutL
		self.NN=np.max(np.abs(self.Ax[:,:,None]-self.Ax[:,None,:]),axis=0)<=1#file format and slicing
		
	def find_Location(self,SAmp):
		#SAmp0=np.clip(SAmp,0,1e12)#raw amplitude estimate--no need (done in detection...)
		###remove 20 percentile
		SAmp-=np.percentile(SAmp,20)#otherwise, I might effectively add stuff...
		#clip to positive values
		SAmp=np.clip(SAmp,0,1e12)/2.#may need to divide by 4 frames (temporal window)
		###Neighboring amplitudes
		#want to look at neighboring amplitudes and have any amplitude smaller than the sum
		#of its neighbors
		Oamp=np.clip(SAmp,0,np.sum(SAmp[:,None]*self.NN,axis=0)/self.OampNorm)
		if self.eCentered:#reverse clipping of central channel
			Oamp[0]=SAmp[0]#want that noisy channels are spatially limited (but would also bias results...)
		###remove 50 percentile (i.e. 37.5 percentile from former baseline)
		ShAmp=np.sum(Oamp)
		CMamp=np.clip(Oamp-np.median(Oamp),0,1e12)
		CM=np.sum(CMamp[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
		#set channels with a distance larger than sqrt(2) to zero/half?
		CMamp2=np.clip(CMamp,0,np.clip(CMamp*(2-np.sqrt(np.sum((CM[:,None]-self.Ax)**2,axis=0))),0,1))
		CM=np.sum(CMamp2[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp2),1e-6,1e12)
		CMamp3=np.clip(CMamp,-0.1*CMamp,CMamp*np.clip((2-np.sqrt(np.sum((CM[:,None]-self.Ax)**2,axis=0))),-0.1,1))
		if np.sum(CMamp3)>np.sum(CMamp)*0.5:
			CMamp=CMamp3
			CM=np.sum(CMamp[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
		else:
			CMamp=np.clip(CMamp,0,CMamp*np.clip((2-np.sqrt(np.sum((CM[:,None]-self.Ax)**2,axis=0))),0,1))
			CM=np.sum(CMamp[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
		ShArea=np.sqrt(np.sum(Oamp[None,:]*(self.Ax-CM[:,None])**2)\
		/np.clip(np.sum(Oamp),1e-6,1e12))
		return CM,CMamp,ShArea
	
	def Iterate(self,z,ts,lSpk):
		###baseline
		if lSpk:
			z0=z[self.IbL]##channels
			z2=z[self.IbL+2]##baseline
			z3=z[self.IbL+3]-z2##amplitudes
			goodCh=(np.abs(z2)<self.Vscale*2000)*self.UseCh[z0]#working channels
		else:
			z0=z[self.Ib]##channels
			z2=z[self.Ib+2]##baseline
			z3=z[self.Ib+3]-z2##amplitudes
			goodCh=(np.abs(z2)<self.Vscale*2000)*self.UseCh[z0]#working channels
		if not all(goodCh):
			if not all(goodCh[:self.Nch0]):
				if lSpk:
					goodCh0=np.arange(self.Nch0,dtype=int)[goodCh[:self.Nch0]]
					zShape=np.reshape(z[:self.LxL],(self.Nch0,self.NcutL0))[goodCh0,4:]
				else:
					goodCh0=np.arange(self.Nch0,dtype=int)[goodCh[:self.Nch0]]
					zShape=np.reshape(z[:self.Lx],(self.Nch0,self.Ncut0))[goodCh0,4:]
				SAmp=np.zeros(self.Nch, dtype=float)
				SAmp[goodCh]=z3[goodCh]*1./np.clip(self.Qdavg[z0[goodCh]],1,1e12)
				###normalize by variance
				zShape*=1./np.clip(self.Qdavg[z0[goodCh0]][:,None],1,1e12)
				for jj in np.nonzero(True-goodCh)[0]:
					SAmp[jj]=np.median(SAmp[self.Mapx[jj]])*self.MapxL[jj]/8.
				CM,CMamp,ShArea=self.find_Location(SAmp)
				#weighted sum of raw traces
				wA=np.sum(zShape*CMamp[goodCh0,None]*1./np.clip(np.sum(CMamp[goodCh0]),1,1e12),axis=0)
				if z0[0]>-1:
					CM+=np.array([z0[0]%64,z0[0]/64])
				else:#find first nonzero
					iind=np.nonzero(z0[0])[0][0]
					CM+=np.array([z0[iind]%64,z0[iind]/64])-self.Ax[iind,:]
			else:
				if lSpk:
					zShape=np.reshape(z[:self.LxL],(self.Nch0,self.NcutL0))[:,4:]
				else:
					zShape=np.reshape(z[:self.Lx],(self.Nch0,self.Ncut0))[:,4:]
				SAmp=np.zeros(self.Nch, dtype=float)
				SAmp[goodCh]=z3[goodCh]*1./np.clip(self.Qdavg[z0[goodCh]],1,1e12)
				###normalize by variance
				zShape*=1./np.clip(self.Qdavg[z0[:self.Nch0]][:,None],1,1e12)
				for jj in np.nonzero(True-goodCh)[0]:
					SAmp[jj]=np.median(SAmp[self.Mapx[jj]])*self.MapxL[jj]/8.
				CM,CMamp,ShArea=self.find_Location(SAmp)
				#weighted sum of raw traces
				wA=np.sum(zShape*CMamp[:self.Nch0,None]*1./np.clip(np.sum(CMamp[:self.Nch0]),1,1e12),axis=0)
				CM+=np.array([z0[0]%64,z0[0]/64])
		else:
			if lSpk:
				zShape=np.reshape(z[:self.LxL],(self.Nch0,self.NcutL0))[:,4:]
			else:
				zShape=np.reshape(z[:self.Lx],(self.Nch0,self.Ncut0))[:,4:]
			SAmp=np.zeros(self.Nch, dtype=float)
			SAmp=z3*1./np.clip(self.Qdavg[z0],1,1e12)
			###normalize by variance
			zShape*=1./np.clip(self.Qdavg[z0[:self.Nch0]][:,None],1,1e12)
			CM,CMamp,ShArea=self.find_Location(SAmp)
			#weighted sum of raw traces
			wA=np.sum(zShape*CMamp[:self.Nch0,None]*1./np.clip(np.sum(CMamp[:self.Nch0]),1,1e12),axis=0)
			CM+=np.array([z0[0]%64,z0[0]/64])
		return CM,CMamp,wA/self.Vscale,ShArea/self.Vscale


### Spike sorting
def readShapesFile(TxtFile, HdfFile, NSpk):
	n=0
	f=h5py.File(HdfFile,'r+')
	Sampling=int(f['Sampling'].value)
	Ncut=int(f['NCut'].value)
	NcutL=int(f['NCutLong'].value)
	PreCut=int(f['PreCut'].value)
	PostCut=int(f['PostCut'].value)
	Reverse=int(f['reverse Detection'].value)
	Ascale=np.abs(int(f['Ascale'].value))
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
	##h=f['GlobalVoltageFluctuations']
	##Bd=np.diff(h['medianVoltage'].value)
	##print Bd.shape
	##Vsbias=h['VFBias'].value
	#invert SortInd
	SortIndInv=np.zeros(SortInd.shape[0],dtype=int)
	SortIndInv[SortInd]=np.arange(SortInd.shape[0])
	##Bsupport=(Sampling*6)/5000+Sampling/1000#need a minimum cutout region
	##BInd=np.zeros((2,Bsupport),dtype=int)#maybe a bit long...
	##BInd[:,:(Sampling)/2000]=1#online baseline estimate
	##BInd[:,Sampling/2000:Sampling/1000]=2+np.arange((Sampling+1000)/2000,dtype=int)[None,:]
	##BInd[0,-(Sampling*6)/5000:]=NcutL+2+np.arange(-(Sampling*6)/5000,0,dtype=int)
	##BInd[1,-(Sampling*6)/5000:]=Ncut+2+np.arange(-(Sampling*6)/5000,0,dtype=int)
	### read shapes
	fName=np.array(['_Shapes.txt','_ShapesX.txt'])
	for iiii in range(2):
		b=file(TxtFile + fName[iiii])
		LF=LocationFinder(iiii==0,Ascale,IgnCh,PreCut,Ncut,NcutL,Qdavg)
		for i in b:
			z=np.array(i.split(),dtype=int)
			if not ((iiii==0)*(int(z[0]) in IgnCh)):
				nInd=SortIndInv[n]
				CM,CMamp,wA,ShArea=LF.Iterate(z,SpkT[nInd],True-PeakL[nInd])
				g['Locations'][nInd,:]=CM[::-1]+0.5#want (y,x) format to be consistent with brainwave
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

class LocationFinderI:
	def __init__(self,eCentered,BInd,Lspike,Ascale,IgnCh,PreCut,Ncut,NcutL,Bd,Qdavg,inserted=False):
		self.Bd=Bd
		self.Qdavg=Qdavg
		self.Vscale=2
		self.PreCut=PreCut
		self.Ncut=Ncut
		self.NcutD=NcutL-self.Ncut
		self.Ascale=Ascale
		self.BInd=BInd
		self.Lspike=Lspike
		self.eCentered=eCentered
		self.UseCh=np.ones(4097,dtype=bool)
		self.UseCh[-1]=False
		self.UseCh[IgnCh]=False
		self.inserted=inserted
		if eCentered:
			self.Mapx=[np.array([5,1,6,2,7,3,8,4]),np.array([6,2,0,4,5]),np.array([1,6,7,3,0]),\
			np.array([4,0,2,7,8]),np.array([5,1,0,3,8]),np.array([1,0,4]),\
			np.array([2,0,1]),np.array([0,2,3]),np.array([4,0,3])]
			self.MapxL=np.array([8,5,5,5,5,3,3,3,3])
			self.Ax=np.array([[0,0,1,0,-1,-1,1,1,-1],[0,-1,0,1,0,-1,-1,1,1]])
			self.Nch=9
			self.OampNorm=np.ones(self.Nch,dtype=float)
			self.OampNorm[1:]+=1.
			self.OampNorm[5:]+=1.
		else:
			self.Mapx=[np.array([4,5,1,2,3,10,11]),np.array([4,5,6,7,2,3,0]),np.array([0,1,3,6,7,8,9]),\
			np.array([0,1,2,8,9,10,11]),np.array([0,1,11,5]),np.array([0,1,4,6]),\
			np.array([1,2,5,7]),np.array([2,1,6,8]),np.array([2,3,7,9]),\
			np.array([2,3,8,10]),np.array([0,3,11,9]),np.array([10,0,3,4])]
			self.MapxL=np.array([7,7,7,7,4,4,4,4,4,4,4,4])
			self.Ax=np.array([[0,1,1,0,0,1,2,2,1,0,-1,-1],[0,0,1,1,-1,-1,0,1,2,2,1,0]])
			self.Nch=12
			self.OampNorm=2*np.ones(self.Nch,dtype=float)
			self.OampNorm[4:]+=1.
		self.NN=np.max(np.abs(self.Ax[:,:,None]-self.Ax[:,None,:]),axis=0)<=1#file format and slicing
		if inserted:
			self.Iweights = np.array([[1000, 190, 250, 250, 250, 250, 190, 190, 190,0],\
			[1000, 208, 250, 300, 250, 214, 174, 208, 174,0],\
			[750, 225, 248, 375, 248, 187, 159, 225, 159,0],\
			[500, 240, 240, 500, 240, 166, 146, 240, 146,0],\
			[1000, 208, 300, 250, 214, 250, 208, 174, 174,0],\
			[1000, 232, 300, 300, 214, 214, 187, 187, 161,0],\
			[750, 258, 297, 375, 213, 187, 168, 198, 149,0],\
			[500, 282, 282, 500, 208, 166, 153, 208, 138,0],\
			[750, 225, 375, 248, 187, 248, 225, 159, 159,0],\
			[750, 258, 375, 297, 187, 213, 198, 168, 149,0],\
			[679, 297, 370, 370, 187, 187, 177, 177, 140,0],\
			[486, 339, 339, 486, 183, 166, 159, 183, 131,0],\
			[500, 240, 500, 240, 166, 240, 240, 146, 146,0],\
			[500, 282, 500, 282, 166, 208, 208, 153, 138,0],\
			[486, 339, 486, 339, 166, 183, 183, 159, 131,0],\
			[414, 414, 414, 414, 163, 163, 163, 163, 123,0]])
			self.Inn = np.argsort(np.array([4, 8, 7, 5, 1, 3, 6, 2, 0]))
			self.Ann = np.array([12, 11, 13, 10, 14, 9, 15, 8, 16, 7, 17, 6, 18, 5, 19])
			self.Iseq = np.array([65, 833, 1601, 2369, 3137, 3905, 577, 1345, 2145, 2913, 3681, 353, 1121, 1889, 2657, 3425, 69, 837, 1605, 2373, 3141, 3909, 581, 1349, 2149, 2917, 3685, 357, 1125, 1893, 2661, 3429, 73, 841, 1609, 2377, 3145, 3913, 585, 1353, 2153, 2921, 3689, 361, 1129, 1897, 2665, 3433, 77, 845, 1613, 2381, 3149, 3917, 589, 1357, 2157, 2925, 3693, 365, 1133, 1901, 2669, 3437, 81, 849, 1617, 2385, 3153, 3921, 593, 1361, 2161, 2929, 3697, 369, 1137, 1905, 2673, 3441, 85, 853, 1621, 2389, 3157, 3925, 597, 1365, 2165, 2933, 3701, 373, 1141, 1909, 2677, 3445, 89, 857, 1625, 2393, 3161, 3929, 601, 1369, 2169, 2937, 3705, 377, 1145, 1913, 2681, 3449, 93, 861, 1629, 2397, 3165, 3933, 605, 1373, 2173, 2941, 3709, 381, 1149, 1917, 2685, 3453, 97, 865, 1633, 2401, 3169, 3937, 609, 1377, 2113, 2881, 3649, 321, 1089, 1857, 2625, 3393, 101, 869, 1637, 2405, 3173, 3941, 613, 1381, 2117, 2885, 3653, 325, 1093, 1861, 2629, 3397, 105, 873, 1641, 2409, 3177, 3945, 617, 1385, 2121, 2889, 3657, 329, 1097, 1865, 2633, 3401, 109, 877, 1645, 2413, 3181, 3949, 621, 1389, 2125, 2893, 3661, 333, 1101, 1869, 2637, 3405, 113, 881, 1649, 2417, 3185, 3953, 625, 1393, 2129, 2897, 3665, 337, 1105, 1873, 2641, 3409, 117, 885, 1653, 2421, 3189, 3957, 629, 1397, 2133, 2901, 3669, 341, 1109, 1877, 2645, 3413, 121, 889, 1657, 2425, 3193, 3961, 633, 1401, 2137, 2905, 3673, 345, 1113, 1881, 2649, 3417, 125, 893, 1661, 2429, 3197, 3965, 637, 1405, 2141, 2909, 3677, 349, 1117, 1885, 2653, 3421])
			self.Template = np.array([[-406, 25101, 32871, 18504, 3709, -4793, -8225, -8875, -8292, -7287, -6233, -5275, -4456, -3777, -3220, -2767, -2397, -2095, -1848, -1643, -1474, -1332, -1214, -1113],\
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
			self.AxCh=self.Ax[0,:]+self.Ax[1,:]*64
			self.Template256=np.zeros((16,256),dtype=int)
			for i in range(16):
				self.Template256[i,np.arange(24)]=self.Template[i,:]
			self.Rxseq=(np.arange(4096)/64/4)*16+((np.arange(4096)/4)%16)
			self.RIseq=np.argsort(self.Rxseq[self.Iseq])
			self.Ix=np.zeros((16,16,self.Nch),dtype=int)
			H=np.zeros((7,7),dtype=int)+9
			H[1:-3,1:-3]=np.reshape(self.Inn,(3,3))
			for ij1 in range(4):
				for ij2 in range(4):
						self.Ix[ij1*4+ij2,:,:]=self.Iweights[:,H[self.Ax[1,:]+1+ij1,self.Ax[0,:]+1+ij2]]
			self.Slist=[]
			for i in range(240):
				self.Slist.append(np.zeros((512,self.Ncut-3,9),dtype=int))
			
	
	def boundedAmplitudes(self,z):
		#bad stuff #ignore channels in list
		return np.all(np.abs(z[:,4:])<self.Vscale*2000,axis=1)*self.UseCh[z[:,0]]#working channels
	
	def ChAmplitudes(self,z,lSpk,Vbias,ts):
		###Amplitudes
		#minimum over Lspike consecutive frames
		#z[:,2:]-=Vsbias[z[:,0],None]*bA2[None,:]
		#bA2=self.Bd[ts-(self.PreCut+1)/2:ts+self.Lspike+self.PreCut][None,:]
		lf=self.Lspike+self.PreCut+self.PreCut/3+4
		li=self.PreCut-self.PreCut/3+4
		if self.Lspike==3:
			#doesn't affect baseline
			Z=np.sum(z[:,li:lf]\
			-self.Bd[ts-(self.PreCut)/3:ts+self.Lspike+self.PreCut/3][None,:]*Vbias[:,None],axis=0)
			lx=np.argmin(2*Z[:-2]+3*Z[1:-1]+2*Z[2:])
			SAmp=-1./7.*(2*z[:,li+lx]+3*z[:,li+lx]+2*z[:,li+lx])
		else:
			Zc=np.cumsum(np.sum(z[:,li-1:lf]-\
			self.Bd[ts-self.PreCut/3-1:ts+self.Lspike+self.PreCut/3][None,:]*Vbias[:,None],axis=0))
			lx=np.argmin(Zc[self.Lspike+1:]-Zc[:-self.Lspike-1]+Zc[self.Lspike:-1]-Zc[1:-self.Lspike],axis=0)
			SAmp=-1./(2.*self.Lspike-2.)*(z[:,li+lx]+z[:,li+lx+self.Lspike-1]\
			+2*np.sum(z[:,li+lx+1:li+lx+self.Lspike-1],axis=1))
		return SAmp
		
	def find_Location(self,SAmp):
		SAmp0=np.clip(SAmp,0,1e12)/2.#raw amplitude estimate
		###remove 20 percentile
		SAmp-=np.percentile(SAmp0,20)#otherwise, I might effectively add stuff...
		#clip to positive values
		SAmp=np.clip(SAmp,0,1e12)#may need to divide by 4 frames (temporal window)
		###Neighboring amplitudes
		#want to look at neighboring amplitudes and have any amplitude smaller than the sum
		#of its neighbors
		Oamp=np.clip(SAmp,0,np.sum(SAmp[:,None]*self.NN,axis=0)/self.OampNorm)
		if self.eCentered:#reverse clipping of central channel
			Oamp[0]=SAmp[0]#want that noisy channels are spatially limited (but would also bias results...)
		###remove 50 percentile (i.e. 37.5 percentile from former baseline)
		ShAmp=np.sum(Oamp)
		CMamp=np.clip(Oamp-np.median(Oamp),0,1e12)
		CM=np.sum(CMamp[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
		#set channels with a distance larger than sqrt(2) to zero/half?
		CMamp2=np.clip(CMamp,0,np.clip(CMamp*(2-np.sqrt(np.sum((CM[:,None]-self.Ax)**2,axis=0))),0,1))
		CM=np.sum(CMamp2[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp2),1e-6,1e12)
		CMamp3=np.clip(CMamp,-0.1*CMamp,CMamp*np.clip((2-np.sqrt(np.sum((CM[:,None]-self.Ax)**2,axis=0))),-0.1,1))
		if np.sum(CMamp3)>np.sum(CMamp)*0.5:
			CMamp=CMamp3
			CM=np.sum(CMamp[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
		else:
			CMamp=np.clip(CMamp,0,CMamp*np.clip((2-np.sqrt(np.sum((CM[:,None]-self.Ax)**2,axis=0))),0,1))
			CM=np.sum(CMamp[None,:]*self.Ax,axis=1)/np.clip(np.sum(CMamp),1e-6,1e12)
		ShArea=np.sqrt(np.sum(Oamp[None,:]*(self.Ax-CM[:,None])**2)\
		/np.clip(np.sum(Oamp),1e-6,1e12))
		return CM,CMamp,ShArea
		
	def Iterate_inserted(self,z,ts,lSpk,Vbias):
		#bA2=self.Bd[ts-self.PreCut:ts-self.PreCut+self.Ncut+self.NcutD*lSpk+1]
		goodCh=self.boundedAmplitudes(z)
		if z[0,0]>-1:
			Iind=((int(z[0,0])/64)%4)*4+(int(z[0,0])%4)
			CMoffset=np.array([z[0,0]%64,z[0,0]/64])
			ICh=self.Rxseq[int(z[0,0])]
		else:#find first nonzero
			iind=np.nonzero(z[:,0])[0][0]
			CMoffset=np.array([z[:,iind]%64,z[:,iind]/64])-self.Ax[iind,:]
			Iind=((int(z[iind])/64-self.Ax[iind,1])%4)*4+(int(z[0]-self.Ax[iind,0])%4)
			ICh=self.Rxseq[int(z[iind,0])-self.AxCh[iind]]
		IT=int(ts-3)
		Itime=(self.RIseq[ICh]-IT)%256#temporal offset (from where to start)
		Atime=((Itime+128)%256-128+IT)%15
		Itrange=np.arange(Itime-2,Itime+self.Ncut-self.PreCut+1)%256
		#print self.Nch,self.Ncut,self.RIseq[ICh]%16, Itrange, Iind, Atime
		z[:self.Nch,self.PreCut+1:self.Ncut+4]-=((self.Template256[self.RIseq[ICh]%16,Itrange])[None,:]\
		*(self.Ix[Iind,self.RIseq[ICh]/16,:])[:,None])/self.Ann[Atime]/500./64.
		if (Iind==5) and (self.Nch==9):
			for ij in range(9):
				if goodCh[ij]:
					self.Slist[self.RIseq[ICh]/16*15+Atime][np.clip(np.array(z[ij,self.PreCut-1:self.Ncut+4]\
					,dtype=int)+400,0,511),np.arange(self.Ncut-self.PreCut+5),ij]+=1
		###baseline
		z[:,1]*=1./self.Ascale
		if lSpk:
			z[:,1:]-=np.median(z[:,self.BInd[0,:]],axis=1)[:,None]
		else:
			z[:,1:]-=np.median(z[:,self.BInd[1,:]],axis=1)[:,None]
		#zChannels[:,noCh]=self.Ax+np.array([z[:,0]%64,z[:,0]/64])[:,None]#no use
		SAmp=np.zeros(self.Nch, dtype=float)
		SAmp[goodCh]=np.array(self.ChAmplitudes(z[goodCh,:],lSpk,Vbias[goodCh],ts),dtype=float)\
		*1./np.clip(self.Qdavg[z[goodCh,0]],1,1e12)
		###normalize by variance
		z[goodCh,4:]/=np.clip(self.Qdavg[z[goodCh,0]][:,None],1,1e12)
		if not all(goodCh):
			for jj in np.nonzero(True-goodCh)[0]:
				z[jj,4:]=np.mean(z[self.Mapx[jj],4:],axis=0)*self.MapxL[jj]/8.
				SAmp[jj]=np.median(SAmp[self.Mapx[jj]])*self.MapxL[jj]/8.
		CM,CMamp,ShArea=self.find_Location(SAmp)
		#weighted sum of raw traces
		wA=np.sum(z[:,4:]*CMamp[:,None]*1./np.clip(np.sum(CMamp),1,1e12),axis=0)
		CM+=CMoffset
		return CM,CMamp,wA,ShArea

### Spike sorting
def readShapesFileInserted(TxtFile, HdfFile, NSpk, Sampling, Lspike=4):
	n=0
	f=h5py.File(HdfFile,'r+')
	Sampling=int(f['Sampling'].value)
	Ncut=int(f['NCut'].value)
	NcutL=int(f['NCutLong'].value)
	PreCut=int(f['PreCut'].value)
	PostCut=int(f['PostCut'].value)
	Reverse=int(f['reverse Detection'].value)
	Ascale=np.abs(int(f['Ascale'].value))
	Bsupport=(Sampling*6)/5000+Sampling/1000#need a minimum cutout region
	BInd=np.zeros((2,Bsupport),dtype=int)#maybe a bit long...
	BInd[:,:(Sampling)/2000]=1#online baseline estimate
	BInd[:,Sampling/2000:Sampling/1000]=4+np.arange((Sampling+1000)/2000,dtype=int)[None,:]
	BInd[0,-(Sampling*6)/5000:]=NcutL+4+np.arange(-(Sampling*6)/5000,0,dtype=int)
	BInd[1,-(Sampling*6)/5000:]=Ncut+4+np.arange(-(Sampling*6)/5000,0,dtype=int)
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
	Bd=np.diff(h['medianVoltage'].value)
	print Bd.shape
	Vsbias=h['VFBias'].value
	#invert SortInd
	SortIndInv=np.zeros(SortInd.shape[0],dtype=int)
	SortIndInv[SortInd]=np.arange(SortInd.shape[0])
	fName=np.array(['_Shapes.txt','_ShapesX.txt'])
	for iiii in range(2):
		b=file(TxtFile + fName[iiii])
		if iiii==0:
			LF=LocationFinderI(iiii==0,BInd,Lspike,Ascale,IgnCh,PreCut,Ncut,NcutL,Bd,Qdavg,inserted=True)
		else:
			LF2=LocationFinderI(iiii==0,BInd,Lspike,Ascale,IgnCh,PreCut,Ncut,NcutL,Bd,Qdavg,inserted=True)
		for i in b:
			z=np.array(i.split(),dtype=int)
			if not ((iiii==0)*(int(z[0]) in IgnCh)):
				nInd=SortIndInv[n]
				z=np.reshape(z,(-1,Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])+4))#know that shape
				if iiii==0:
					CM,CMamp,wA,ShArea=LF.Iterate_inserted(z,SpkT[nInd],True-PeakL[nInd],Vsbias[z[:,0]])
				else:
					CM,CMamp,wA,ShArea=LF2.Iterate_inserted(z,SpkT[nInd],True-PeakL[nInd],Vsbias[z[:,0]])
				g['Locations'][nInd,:]=CM[::-1]+0.5#want (y,x) format to be consistent with brainwave
				g['ShAmpX'][nInd]=np.sum(CMamp)/4.#might be less noisy... and will subtract a larger baseline for wide spikes
				g['ShArea'][nInd]=ShArea/4.
				g['Shapes'][nInd,:Ncut*PeakL[nInd]+NcutL*(1-PeakL[nInd])]=wA/4.
				n+=1
		b.close()
	MedianShape=np.zeros((15,16,Ncut-3,9))
	MedianShapeDev=np.zeros((15,16,Ncut-3,9))
	Q=((np.arange(Ncut-3,dtype=int)-2)%256)
	Template256=LF.Template256
	Slist=LF.Slist
	TemplateMean=np.mean(Template256[:,Q],axis=0)
	for i in range(240):
		a=np.sum(Slist[i][:,0,:],axis=0)
		for j in range(9):
			MedianShape[i%15,i/15,:,j]=(np.argmax(np.diff(np.concatenate((np.zeros((1,Ncut-3))\
			,1*(np.cumsum(Slist[i][:,:,j],axis=0)>a[j]/2)),axis=0),axis=0),axis=0)-400)/4.
			MedianShapeDev[i%15,i/15,:,j]=-MedianShape[i%15,i/15,:,j]-(TemplateMean\
			*LF.Ix[5,i/15,j])*1./LF.Ann[i%15]/1000./LF.Ascale
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
	AmplitudesB=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/2.
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
		Ispikes=RepolarizingSpikes0.copy()
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

def IsolatedSpikes(HdfFile, IncludeLongSpikes=True, DFrames=2, MaxDist=1.):
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
	Amplitudes=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/2.#+1e-6*scipy.rand(Times.shape())#not like that!
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
	Amplitudes=ShAmpX*np.clip((1.-np.abs(ShArea-1.)),0,1)/2.
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
