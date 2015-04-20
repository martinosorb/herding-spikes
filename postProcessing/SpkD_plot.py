import numpy as np
import pylab
import os.path
import scipy.stats
import matplotlib as mplt
import h5py
import scipy.cluster
#import mpl_toolkits.axisartist as AA

font_size = 10
mplt.rcParams['axes.titlesize'] = font_size+2
mplt.rcParams['xtick.labelsize'] = font_size
mplt.rcParams['ytick.labelsize'] = font_size
mplt.rcParams['axes.labelsize'] = font_size+1
mplt.rcParams['legend.fontsize'] = font_size-2
mplt.rcParams['font.size'] = font_size+2
mplt.rcParams['font.family']='sans-serif'
mplt.rcParams['font.sans-serif']=['Arial']
mplt.rcParams['image.interpolation']='none'
mplt.rcParams['ps.fonttype']=42
mplt.rcParams['savefig.pad_inches']=0.
mplt.rcParams['savefig.bbox']='tight'
mplt.rcParams['savefig.dpi']=300

def Rasterplot(HdfFile, FileName, FieldName=''\
, tMin=0, tMax=20, yMin=0, yMax=64, cMin=0, cMax=1., cRes=20, cLabel=''\
, RepolSpikes=True, LongSpikes=True, axis=0, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.5, onlineVersion=False):
	assert cMin<cMax
	assert tMin<tMax
	assert cRes>=2
	A=np.array(['y-position','x-position'])
	g=h5py.File(HdfFile,'r')
	fig=pylab.figure(1,figsize=(4,3))
	ax=fig.add_axes([0.13, 0.13, 0.7, 0.8],axisbg=bgColor)
	if LongSpikes and RepolSpikes and g['IncludeLongSpikes'].value:
		ax0=fig.add_axes([0.84, 0.13, 0.015, 0.8])
		ax1=fig.add_axes([0.855, 0.13, 0.015, 0.8])
	else:
		ax1=fig.add_axes([0.845, 0.13, 0.025, 0.8])
	ax.tick_params(pad=4)
	Sampling=g['Sampling'].value
	Times=g['Times'].value
	tInd=(Times>tMin*Sampling)*(Times<tMax*Sampling)
	Times=Times[tInd]
	Loc=g['Locations'].value[tInd,axis]
	if 'CorrelationAnalysis/Probability' in g:
		p=g['CorrelationAnalysis/Probability'].value[tInd]
	else:
		p=1.
	Amp=g['Amplitudes'].value[tInd]
	pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	Times=Times[pAmpInd]
	Loc=Loc[pAmpInd]
	if not (FieldName==''):#plot the event property
		cQ=np.clip(np.array((g[FieldName].value[tInd][pAmpInd]-cMin)*cRes*1./(cMax-cMin),dtype=int),0,cRes-1)
	else:#make a 2d plot
		cQ=np.clip(np.array((g['Locations'].value[tInd,1-axis][pAmpInd]-cMin)*cRes*1./(cMax-cMin),dtype=int),0,cRes-1)
	if g['IncludeLongSpikes'].value:
		sInd=g['RepolarizingSpikes'].value[tInd][pAmpInd]
		if not RepolSpikes:
			sInd=True-sInd
			LongSpikes=False
		if LongSpikes:
			Times0=Times[True-sInd]
			Loc0=Loc[True-sInd]
			cQ0=cQ[True-sInd]
			for i in range(cRes):
				ax.plot(Times0[cQ0==i]*1./Sampling,Loc0[cQ0==i],',',color=LongSpikesCmap(i*1./(cRes-1)), alpha=alpha)
		Times=Times[sInd]
		Loc=Loc[sInd]
		cQ=cQ[sInd]
		for i in range(cRes):
			ax.plot(Times[cQ==i]*1./Sampling,Loc[cQ==i],',',color=SpikesCmap(i*1./(cRes-1)), alpha=alpha)
		if LongSpikes:
			cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=LongSpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
			cbar0.set_ticks(np.array([]))
	elif onlineVersion:
		Loc2=g['Locations'].value[tInd,1-axis][pAmpInd]
		for i in range(cRes):
			ax.plot(Times[cQ==i]*1./Sampling,Loc[cQ==i]+Loc2[cQ==i]/64.,',',color=SpikesCmap(i*1./(cRes-1)), alpha=alpha)	
	else:
		for i in range(cRes):
			ax.plot(Times[cQ==i]*1./Sampling,Loc[cQ==i],',',color=SpikesCmap(i*1./(cRes-1)), alpha=alpha)	
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	if not (cLabel==''):
		cbar1.set_label(cLabel,fontsize='small')
	elif not (FieldName==''):
		cbar1.set_label(FieldName,fontsize='small')
	else:
		cbar1.set_label(A[1-axis],fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	ax.set_xlabel('time/s')
	ax.set_ylabel(A[axis])
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(tMin,tMax)
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	print Times.shape
	return
	
def Scatterplot(HdfFile, FileName, FieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=10, cRes=20, cLabel=''\
, RepolSpikes=True, LongSpikes=True, ampThreshold=0, pThreshold=0, bgColor=pylab.cm.gray(0.5)\
, SpikesCmap=pylab.cm.hot_r, LongSpikesCmap=pylab.cm.cool, alpha=0.5):
	assert cMin<cMax
	assert xMin<xMax
	assert cRes>=2
	A=np.array(['y-position','x-position'])
	g=h5py.File(HdfFile,'r')
	fig=pylab.figure(1,figsize=(4,3.5))
	ax=fig.add_axes([0.13, 0.13, 0.7, 0.8],axisbg=bgColor)
	if LongSpikes and RepolSpikes and g['IncludeLongSpikes'].value:
		ax0=fig.add_axes([0.84, 0.13, 0.015, 0.8])
		ax1=fig.add_axes([0.855, 0.13, 0.015, 0.8])
	else:
		ax1=fig.add_axes([0.845, 0.13, 0.025, 0.8])
	ax.tick_params(pad=4)
	Sampling=g['Sampling'].value
	if 'CorrelationAnalysis/Probability' in g:
		p=g['CorrelationAnalysis/Probability'].value
	else:
		p=1.
	Amp=g['Amplitudes'].value
	pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	Loc=g['Locations'].value[pAmpInd,:]
	cQ=np.clip(np.array((g[FieldName].value[pAmpInd]-cMin)*cRes*1./(cMax-cMin),dtype=int),0,cRes-1)
	if g['IncludeLongSpikes'].value:
		sInd=g['RepolarizingSpikes'].value[pAmpInd]
		if not RepolSpikes:
			sInd=True-sInd
			LongSpikes=False
		if LongSpikes:
			Loc0=Loc[True-sInd,:]
			cQ0=cQ[True-sInd]
			for i in range(cRes):
				ax.plot(Loc0[cQ0==i,1],Loc0[cQ0==i,0],',',color=LongSpikesCmap(i*1./(cRes-1)), alpha=alpha)
		Loc=Loc[sInd,:]
		cQ=cQ[sInd]
		for i in range(cRes):
			ax.plot(Loc[cQ==i,1],Loc[cQ==i,0],',',color=SpikesCmap(i*1./(cRes-1)), alpha=alpha)
		if LongSpikes:
			cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=LongSpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
			cbar0.set_ticks(np.array([]))
	else:
		for i in range(cRes):
			ax.plot(Loc[cQ==i,1],Loc[cQ==i,0],',',color=SpikesCmap(i*1./(cRes-1)), alpha=alpha)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	if not (cLabel==''):
		cbar1.set_label(cLabel,fontsize='small')
	elif not (FieldName==''):
		cbar1.set_label(FieldName,fontsize='small')
	else:
		cbar1.set_label(A[1-axis],fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#ax.set_xlabel('time/s')
	#ax.set_ylabel(A[axis])
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	print Loc.shape
	return

def Densityplot(HdfFile, FileName, FieldName='', ProbFieldName='CorrelationAnalysis/Probability'\
, LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=1, cMax=10000, Res=8, logScale=True, cLabel=''\
, RepolSpikes=True, LongSpikes=True, ampThreshold=0, pThreshold=0, cThreshold=-1e12\
, nThreshold=0, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r, alpha=1.):
	assert cMin<cMax
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(4,3.5))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=4)
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	if not (ProbFieldName==''):
		p=g[ProbFieldName].value
	else:
		p=1.
	Amp=g[AmpFieldName].value
	cQb=False
	if not (FieldName==''):#plot the event property
		cQ=g[FieldName].value
		cQb=True
		pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)*(cQ>=cThreshold)
		cQ=cQ[pAmpInd]
	else:
		pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	Loc=g[LocFieldName].value[pAmpInd,:]
	if g['IncludeLongSpikes'].value:
		sInd=g['RepolarizingSpikes'].value[pAmpInd]
		if not RepolSpikes:
			sInd=True-sInd
			LongSpikes=False
		if LongSpikes and RepolSpikes:
			H=np.histogram2d(Loc[:,0],Loc[:,1]\
			,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
			H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
			if cQb:
				Hq=np.histogram2d(Loc[:,0],Loc[:,1]\
				,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res)\
				,weights=cQ)[0]
				Hq=Hq[1:,1:]+Hq[:-1,1:]+Hq[:-1,:-1]+Hq[1:,:-1]
		else:
			H=np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
			,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
			H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
			if cQb:
				Hq=np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
				,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res)\
				,weights=cQ[sInd])[0]
				Hq=Hq[1:,1:]+Hq[:-1,1:]+Hq[:-1,:-1]+Hq[1:,:-1]
	else:
		H=np.histogram2d(Loc[:,0],Loc[:,1]\
		,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
		H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
		if cQb:
			Hq=np.histogram2d(Loc[:,0],Loc[:,1]\
			,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res)\
			,weights=cQ)[0]
			Hq=Hq[1:,1:]+Hq[:-1,1:]+Hq[:-1,:-1]+Hq[1:,:-1]
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	if not (cLabel==''):
		cbar1.set_label(cLabel,fontsize='small')
	elif not (FieldName==''):
		cbar1.set_label(FieldName,fontsize='small')
	else:
		cbar1.set_label('Event count',fontsize='small')
	if logScale and not cQb:
		cbar1.set_ticks((np.arange(int(np.ceil(np.log10(cMin))),int(np.log10(cMax))+1)-np.log10(cMin))/np.log10(cMax*1./cMin))
		cbar1.set_ticklabels(10**np.arange(int(np.ceil(np.log10(cMin))),int(np.log10(cMax))+1))
	else:
		cbar1.set_ticks(np.array([0,0.5,1.0]))
		cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	if cQb:
		kMax=int(np.ceil(1./alpha))
		if (kMax>1):
			for k in range(kMax):
				ax.imshow(np.ma.array(Hq*1./np.clip(H,1,1e12),mask=True-((H>=(2**(kMax-k-1)))*(H>nThreshold)))\
				,cmap=SpikesCmap, aspect='equal'\
				,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax), alpha=alpha)
				ax0=fig.add_axes([0.845, 0.13+(0.025*k)/kMax, 0.025/kMax, 0.8])
				ax0.set_frame_on(False)
				ax0.set_yticks(np.array([]))
				ax0.set_xticks(np.array([]))
				ax0.imshow(np.outer(np.ones(100),np.ones(5)),cmap=pylab.cm.gray,origin='lower',aspect='auto')
				ax0.imshow(np.outer(np.linspace(0,1,100),np.ones(5)),alpha=alpha\
				,cmap=LongSpikesCmap,origin='lower',aspect='auto')
		else:
			ax.imshow(np.ma.array(Hq*1./np.clip(H,1,1e12),mask=True-((H>=1)*(H>nThreshold)))\
			,vmin=cMin, vmax=cMax,cmap=SpikesCmap, aspect='equal'\
			,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax), alpha=alpha)
	elif logScale:
		ax.imshow(np.ma.array(H*Res**2/4.,mask=True-((H>=1)*(H>nThreshold))),norm=mplt.colors.LogNorm(cMin,cMax)\
		,cmap=SpikesCmap, aspect='equal'\
		,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	else:
		ax.imshow(H*Res**2/4.,vmin=cMin, vmax=cMax, cmap=SpikesCmap, aspect='equal'\
		,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	print Loc.shape
	return


def Matrixplot(HdfFile, FileName, FieldName='CorrelationAnalysis/Noise'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., Res=3, logScale=False, cLabel=''\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	assert cMin<cMax
	assert xMin<xMax
	fig=pylab.figure(1,figsize=(4,3.5))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=4)
	g=h5py.File(HdfFile,'r')
	cQ=g[FieldName].value
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	if not (cLabel==''):
		cbar1.set_label(cLabel,fontsize='small')
	else:
		cbar1.set_label(FieldName,fontsize='small')
	if logScale:
		cbar1.set_ticks((np.arange(int(np.ceil(np.log10(cMin))),int(np.log10(cMax))+1)-np.log10(cMin))/np.log10(cMax*1./cMin))
		cbar1.set_ticklabels(10**np.arange(int(np.ceil(np.log10(cMin))),int(np.log10(cMax))+1))
		ax.imshow(np.reshape(cQ[:,0],(64*Res,64*Res))[yMin*Res:yMax*Res,xMin*Res:xMax*Res],norm=mplt.colors.LogNorm(cMin,cMax)\
		,cmap=SpikesCmap, aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	else:
		cbar1.set_ticks(np.array([0,0.5,1.0]))
		cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
		ax.imshow(np.reshape(cQ[:,0],(64*Res,64*Res))[yMin*Res:yMax*Res,xMin*Res:xMax*Res],vmin=cMin, vmax=cMax\
		,cmap=SpikesCmap, aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return

def	Clusterplot(HdfFile, FileName\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=10, cMax=10000, Res=12\
, RepolSpikes=True, LongSpikes=True, ampThreshold=0, pThreshold=0, bgColor='w'\
, ClusterCmap=pylab.cm.hsv):
	fig=pylab.figure(1,figsize=(4,3.5))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=4)
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	if 'CorrelationAnalysis/Probability' in g:
		p=g['CorrelationAnalysis/Probability'].value
	else:
		p=1.
	Amp=g['Amplitudes'].value
	pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	Loc=g['Locations'].value[pAmpInd,:]
	if g['IncludeLongSpikes'].value:
		sInd=g['RepolarizingSpikes'].value[pAmpInd]
		if not RepolSpikes:
			sInd=True-sInd
			LongSpikes=False
		if LongSpikes and RepolSpikes:
			H=np.histogram2d(Loc[:,0],Loc[:,1]\
			,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
			H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
		else:
			H=np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
			,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
			H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
	else:
		H=np.histogram2d(Loc[:,0],Loc[:,1]\
		,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
		H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
	nBins=g['Cluster/Parameter/nBins'].value
	Xmap=g['Cluster/CBoundaries'].value.flatten()
	Xnew=g['Cluster/CAreaMatrix'].value.flatten()
	NChannels=g['Cluster/NCount'].value.shape[0]
	Yi=np.argsort(scipy.rand(NChannels))
	Y=Yi[Xnew]
	ClusterCmap.set_bad(bgColor,alpha=0)
	ax.imshow(np.ma.array(H*Res**2/4.,mask=H<1),norm=mplt.colors.LogNorm(cMin,cMax)\
	,cmap=pylab.cm.gray_r, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	CImg=np.reshape(np.ma.array(Y,mask=((Xmap==3)+(Xmap==0)+(Xnew==0))),(nBins,nBins))
	ax.imshow(CImg[(xMin*nBins)/64:(xMax*nBins)/64,:][:,(yMin*nBins)/64:(yMax*nBins)/64]\
	,cmap=ClusterCmap, aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	cbar =mplt.colorbar.ColorbarBase(ax=ax1,cmap=ClusterCmap,norm=mplt.colors.Normalize(vmin=0,vmax=NChannels))
	cbar.set_label(r'cluster id')
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return
	
def Shapesplot(HdfFile, FileName\
, yMin=-20, yMax=10, yBins=100, xVal=23.5, yVal=23.5, cMin=1, cMax=1000, radius=0.5, nMax=10000\
, RepolSpikes=True, LongSpikes=True, ampThreshold=0, pThreshold=0, bgColor='w'\
, SpikesCmap=pylab.cm.hot_r):
	fig=pylab.figure(1,figsize=(4,3.5))
	ax=fig.add_axes([0.15, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.865, 0.13, 0.025, 0.8])
	ax.tick_params(pad=4)
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	PeakIndex=g['PeakIndex'].value
	if 'CorrelationAnalysis/Probability' in g:
		p=g['CorrelationAnalysis/Probability'].value
	else:
		p=1.
	Amp=g['Amplitudes'].value
	Shapes=g['Shapes'].value
	pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	Loc=g['Locations'].value[pAmpInd,:]
	if g['IncludeLongSpikes'].value:
		sInd=g['RepolarizingSpikes'].value[pAmpInd]
		if not RepolSpikes:
			sInd=True-sInd
			LongSpikes=False
		if LongSpikes and RepolSpikes:
			Ind =np.nonzero(np.sum((Loc-np.array([yVal,xVal])[None,:])**2,axis=1)<radius**2)[0]
		else:
			Ind =np.nonzero((np.sum((Loc-np.array([yVal,xVal])[None,:])**2,axis=1)<radius**2)*sInd)[0]
	else:
		Ind =np.nonzero(np.sum((Loc-np.array([yVal,xVal])[None,:])**2,axis=1)<radius**2)[0]
	if len(Ind)>nMax:
		Ind=Ind[np.argsort(np.argsort(np.sum((Loc[Ind,:]-np.array([yVal,xVal])[None,:])**2,axis=1)))<nMax]
	H=np.histogram2d(Shapes[Ind,:].flatten(),np.tile(np.arange(Shapes.shape[1]),len(Ind))\
	,bins=(np.linspace(yMin,yMax,yBins+1),np.arange(Shapes.shape[1]+1)))[0]
	SpikesCmap.set_bad(bgColor)
	ax.imshow(np.ma.array(H,mask=H<1),norm=mplt.colors.LogNorm(cMin,cMax)\
	,cmap=SpikesCmap, aspect='auto',interpolation='none',origin='lower'\
	,extent=(-(PeakIndex+0.5)*1000./Sampling,(Shapes.shape[1]-PeakIndex-0.5)*1000./Sampling,yMin,yMax))
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_ticks((np.arange(int(np.ceil(np.log10(cMin))),int(np.log10(cMax))+1)-np.log10(cMin))/np.log10(cMax*1./cMin))
	cbar1.set_ticklabels(10**np.arange(int(np.ceil(np.log10(cMin))),int(np.log10(cMax))+1))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(-(PeakIndex+0.5)*1000./Sampling,(Shapes.shape[1]-PeakIndex-0.5)*1000./Sampling)
	ax.set_xlabel('time/ms')
	ax.set_ylabel('avg. ADCcounts')
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return

def GlobalVoltageplot(HdfFile, FileName, cMin=-0.2, cMax=0.3, CorrCmap=pylab.cm.hsv):
	fig=pylab.figure(1,figsize=(8,4))
	ax=fig.add_axes([0.1, 0.45, 0.35, 0.5])
	ax1=fig.add_axes([0.88, 0.4, 0.015, 0.55])
	#ax2=fig.add_axes([0.0625, 0.1, 0.125, 0.25])
	#ax2=fig.add_axes([0.5, 0.4, 0.125, 0.25])
	#ax2=fig.add_axes([0.5, 0.7, 0.125, 0.25])
	xPos=np.array([0.5,0.625,0.75,0.5,0.625,0.75,0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.8125])
	yPos=np.array([0.7,0.7,0.7,0.4,0.4,0.4,0.1,0.1,0.1,0.1,0.1,0.1,0.1])
	ax.tick_params(pad=4)
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	recordedChannels=g['recordedChannels'].value
	B=g['GlobalVoltageFluctuations/medianVoltage'].value
	sumSqVglobal=g['GlobalVoltageFluctuations/sumSqVglobal'].value
	sumSqVchannels=g['GlobalVoltageFluctuations/sumSqVchannels'].value
	sumVproduct=g['GlobalVoltageFluctuations/sumVproduct'].value
	sumVFBias=g['GlobalVoltageFluctuations/VFBias'].value
	sumSqVFBias=g['GlobalVoltageFluctuations/StdVFBias'].value
	Cavg=np.zeros((13,5))
	p=np.array([10,30,50,70,90])
	#plot correlations
	for i in range(13):
		C=np.zeros(4096)
		C[recordedChannels]=sumVproduct[:,i]*1./np.clip(np.sqrt(sumSqVchannels*sumSqVglobal),1,None)
		for j in range(5):
			Cavg[i,j]=np.percentile(C[recordedChannels],p[j])
		ax2=fig.add_axes([xPos[i], yPos[i], 0.125, 0.25])
		ax2.imshow(np.reshape(C,(64,64)),cmap=CorrCmap,aspect='equal'\
		,vmin=cMin,vmax=cMax,interpolation='none',origin='lower',extent=(0,64,0,64))
		ax2.set_xticks(np.array([]))
		ax2.set_yticks(np.array([]))
		ax2.set_xlabel('%4.2f ms'%((i-6)*1000./Sampling))
		cbar =mplt.colorbar.ColorbarBase(ax=ax1,cmap=CorrCmap,norm=mplt.colors.Normalize(vmin=cMin,vmax=cMax))
		cbar.set_label(r'correlation coefficient')
	k=np.array(['k:','k--','k-','k--','k:'])
	for j in range(5):
		ax.plot(np.arange(-6,7)*1000./Sampling,Cavg[:,j],k[j],label='%i'%p[j]+r' %')
	ax.set_xlabel('time lag/ms')
	ax.set_ylabel('correlation coefficient')
	ax.legend()
	ax.set_ylim(cMin,cMax)
	ax.set_xlim(-6000./Sampling,6000./Sampling)
	fig.text(0.1,0.9,'s0=%4.1f muV'%(np.sqrt(sumSqVglobal*1./Sampling/tMax)*4.125/2.048))
	fig.text(0.1,0.82,'s=%4.1f muV'%(np.sqrt(np.mean(sumSqVchannels)*1./Sampling/tMax)*4.125/2.048))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return


def PrintStats(HdfFile,StatsFile, onlineVersion=False):
	if onlineVersion:
		HighThr=8
	else:
		HighThr=2.5
	g=h5py.File(HdfFile,'r')
	f=file(StatsFile,'w')
	f.write('Sampling rate: %i \n'%(g['Sampling'].value))
	f.write('Duration of recording: %i s \n'%(g['tMax'].value))
	f.write('Number of recorded channels: %i \n \n'%(len(g['recordedChannels'].value)))
	if not onlineVersion:
		f.write('Amplitude thresholds: %i %i \n'%(g['AmplitudeThresholds'].value[0], g['AmplitudeThresholds'].value[1]))
		f.write('Repolarization threshold: %i \n'%(g['RepolarizationThreshold'].value))
		f.write('Length of cutouts: %i frames \n'%(g['lenCutouts'].value))
		f.write('Avg. variability estimate of channels: %4.2f \n'%(np.mean(g['ChannelVariability'].value)))
	f.write('Detection \n')
	if 'RawEvents' in g:
		f.write('Number of detected events: %i \n'%(len(g['RawEvents/PreSelectedEvents'].value)))
		f.write('-   repolarizing events: %i \n'%(np.sum(g['RawEvents/RepolarizingSpikes'].value)))
		f.write('-   after discarding noisy channels: %i \n'%(len(g['RawEvents/IsolatedSpikes'].value)))
		f.write('-   after removing multiple detections of the same spike: %i \n '%(np.sum(g['RawEvents/IsolatedSpikes'].value)))
		f.write('-   fraction of unique events: %4.3f \n'%(np.mean(1.*g['RawEvents/IsolatedSpikes'].value)))
	f.write('-   Repolarizing events \n')
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	f.write('-   after removing multiple detections: %i \n '%(np.sum(SingleSpk)))
	f.write('-   fraction of repolarizing events: %4.3f \n \n'%(np.mean(1.*SingleSpk)))
	f.write('-   high amplitude repolarizing events (%2.1f): %i \n \n'\
	%(HighThr,np.sum((g['Amplitudes'].value>=HighThr)*g['RepolarizingSpikes'].value)))
	
	f.write('Correlation Analysis \n')
	if not onlineVersion:
		NChannels=g['CorrelationAnalysis/Parameter/Resolution'].value**2
	else:
		NChannels=64**2
	NSpikes=np.histogram2d(g['Units'].value,(True-SingleSpk)\
	,bins=(np.arange(NChannels+1),np.arange(3)))[0]
	f.write('Number of noisy spikes %i \n'%(np.sum(g['CorrelationAnalysis/Noise'].value*NSpikes)))
	f.write('Number of true spikes %i \n'%(np.sum((1.-g['CorrelationAnalysis/Noise'].value)*NSpikes)))
	f.write('Fraction of noisy spikes %4.3f \n \n'%(np.sum(g['CorrelationAnalysis/Noise'].value*NSpikes)\
	*1./np.sum(NSpikes)))
	if not onlineVersion:
		f.write('Number of noisy spikes (repolarizing) %i \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,0])))
		f.write('Number of true spikes (repolarizing) %i \n'%(np.sum(((1.-g['CorrelationAnalysis/Noise'].value)*NSpikes)[:,0])))
		f.write('Fraction of noisy spikes (repolarizing) %4.3f \n \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,0])\
		*1./np.sum(NSpikes[:,0])))
		f.write('Number of noisy long spikes %i \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,1])))
		f.write('Fraction of noisy long spikes %4.3f \n \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,1])\
		*1./np.sum(NSpikes[:,1])))
		f.write('Cluster \n')
		f.write('Number of clusters %i \n'%(len(g['Cluster/NCount'].value)-1))
		f.write('Number of events in clusters %i \n'%(np.sum(g['Cluster/NCount'].value[1:])))
	f.close()
	g.close()
	return

def PeriEventActivityPlot(HdfFile, FileName):
	g=h5py.File(HdfFile,'r+')
	f=g['PeriEventActivity']
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	Na=int(f['Na'].value)
	Ns=int(f['Ns'].value)
	Res=int(f['Res'].value)
	Aavg=f['Aavg'].value
	Aloc1=f['Aloc'].value
	Davg=f['Davg'].value
	Dloc1=f['Dloc'].value
	Nloc1=f['Nloc'].value
	nBins=Res*64
	#convert to colour
	a=np.array([[1,0,0],[0.5,0.5,0],[0,1,0],[0,0.5,0.5],[0,0,1]])
	cDloc1=np.dot(Dloc1,a)
	cDloc1=cDloc1*1./np.clip(np.sum(cDloc1,axis=2),1e-6,1e12)[:,:,None]
	#need to compute something related to significance
	cDavg=np.clip(np.dot((Dloc1-Davg)*np.sqrt(Nloc1)[:,:,None]*Ns/3.,a)+0.5,0,1)
	b=np.array([[1,0,0],[0.5,0.5,0],[0,1,0],[0,0.5,0.5],[0,0,1],[0.5,0,0.5]])
	cAloc1=np.dot(Aloc1,b)
	cAloc1=cAloc1*1./np.clip(np.sum(cAloc1,axis=2),1e-6,1e12)[:,:,None]
	cAavg=np.clip(np.dot((Aloc1-Aavg)*np.sqrt(Nloc1)[:,:,None]*Ns/3.,b)+0.5,0,1)
	#colorbar
	A=np.angle(np.linspace(-1,1,200)[None,:]*np.ones(200)[:,None]\
	+1j*np.linspace(-1,1,200)[:,None]*np.ones(200)[None,:])*3./np.pi-0.5
	Ac=np.zeros((200,200,6))
	for i in range(6):
		Ac[:,:,i]+=(1.-A%1)*(np.array(A%6,dtype=int)==i)
		Ac[:,:,i]+=(A%1)*(np.array((A+1)%6,dtype=int)==i)
	C=np.abs(np.linspace(-1,1,200)[None,:]*np.ones(200)[:,None]\
	+1j*np.linspace(-1,1,200)[:,None]*np.ones(200)[None,:])
	Cbar=(C[:,:,None]*np.dot(Ac,b)+1./3.*(1.-C[:,:,None]))*((C<1.)[:,:,None])
	#plotting----------------------------------------------------------------------------
	fig2=pylab.figure(2,figsize=(6.4,3.2))
	fig3=pylab.figure(3,figsize=(6.5,3))
	#pylab.subplots_adjust(left=0.12, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3)
	ax0=fig2.add_axes([0.05, 0.1, 0.38, 0.76])
	ax3=fig3.add_axes([0.05, 0.08, 0.35, 0.78])
	#ax2=fig3.add_axes([0.63, 0.54, 0.25, 0.38])
	ax4=fig3.add_axes([0.45, 0.08, 0.35, 0.78])
	ax1=fig2.add_axes([0.48, 0.1, 0.38, 0.76])
	ax5=fig3.add_axes([0.81, 0.08, 0.18, 0.45])
	ax6=fig2.add_axes([0.87, 0.1, 0.02, 0.76])
	#ax7=fig3.add_axes([0.92, 0.05, 0.015, 0.38])
	ax0.tick_params(pad=4)
	ax1.tick_params(pad=4)
	#ax2 .tick_params(pad=4)
	ax3.tick_params(pad=4)
	ax4.tick_params(pad=4)
	#ax5.tick_params(pad=4)
	fig2.text(0.01,0.91,'A',fontsize='xx-large')
	fig2.text(0.41,0.91,'B',fontsize='xx-large')
	fig3.text(0.01,0.91,'A',fontsize='xx-large')
	fig3.text(0.41,0.91,'B',fontsize='xx-large')
	#fig3.text(0.31,0.45,'E',fontsize='xx-large')
	#fig3.text(0.61,0.45,'F',fontsize='xx-large')
	ax0.imshow(cDloc1,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,nBins-Ns+1,0,nBins-Ns+1))
	ax0.set_xlim(-2,nBins-2)
	ax0.set_ylim(-2,nBins-2)
	ax0.set_xticks(np.array([0,nBins/2,nBins])-2)
	ax0.set_xticklabels(np.array([0,32,64]))
	ax0.set_yticks(np.array([0,nBins/2,nBins])-2)
	ax0.set_yticklabels(np.array([0,32,64]))
	ax3.imshow(cAloc1,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,nBins-Ns+1,0,nBins-Ns+1))
	ax3.set_xlim(-2,nBins-2)
	ax3.set_ylim(-2,nBins-2)
	ax3.set_xticks(np.array([0,nBins/2,nBins])-2)
	ax3.set_xticklabels(np.array([0,32,64]))
	ax3.set_yticks(np.array([0,nBins/2,nBins])-2)
	ax3.set_yticklabels(np.array([0,32,64]))
	ax1.imshow(cDavg,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,nBins-Ns+1,0,nBins-Ns+1))
	ax1.set_xlim(-2,nBins-2)
	ax1.set_ylim(-2,nBins-2)
	ax1.set_xticks(np.array([0,nBins/2,nBins])-2)
	ax1.set_xticklabels(np.array([0,32,64]))
	ax1.set_yticks(np.array([0,nBins/2,nBins])-2)
	ax1.set_yticklabels(np.array([0,32,64]))
	ax4.imshow(cAavg,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,nBins-Ns+1,0,nBins-Ns+1))
	ax4.set_xlim(-2,nBins-2)
	ax4.set_ylim(-2,nBins-2)
	ax4.set_xticks(np.array([0,nBins/2,nBins])-2)
	ax4.set_xticklabels(np.array([0,32,64]))
	ax4.set_yticks(np.array([0,nBins/2,nBins])-2)
	ax4.set_yticklabels(np.array([0,32,64]))
	rgbadict = {'red': ((0.0, 0.75, 0.75),
                 (1./12., 1.0, 1.0),
                 (5./12., 0.0, 0.0),
                 (9./12., 0.0, 0.0),
                 (1.0, 0.75, 0.75)),
         'green': ((0.0, 0.0, 0.0),
                 (1./12., 0.0, 0.0),
                 (5./12., 1.0, 1.0),
                 (9./12., 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.25, 0.25),
                 (1./12., 0.0, 0.0),
                 (5./12., 0.0, 0.0),
                 (9./12., 1.0, 1.0),
                 (1.0, 0.25, 0.25))}
	rgbacmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbadict,256)
	rgbdict = {'red': ((0.0, 1.0, 1.0),
                 (0.5, 0.0, 0.0),
                 (1.0, 0.0, 0.0)),
         'green': ((0.0, 0.0, 0.0),
                 (0.5, 1.0, 1.0),
                 (1.0, 0.0, 0.0)),
         'blue': ((0.0, 0.0, 0.0),
                 (0.5, 0.0, 0.0),
                 (1.0, 1.0, 1.0))}
	rgbcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbdict,256)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax6,cmap=rgbcmap,norm=mplt.colors.Normalize(2*0.042,26*0.042))
	cbar0.set_label(r'avg. distance/mm')
	ax5.imshow(Cbar,aspect='equal'\
	,interpolation='none',origin='lower')
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	#cbar0 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=rgbacmap,norm=mplt.colors.Normalize(0,2.*np.pi))
	#cbar0.set_label(r'avg. angle',fontsize='x-small')
	pylab.figure(2)
	pylab.savefig(FileName+'_dist.png')
	pylab.savefig(FileName+'_dist.eps')
	pylab.figure(3)
	pylab.savefig(FileName+'_angle.png')
	pylab.savefig(FileName+'_angle.eps')
	pylab.close()
	return
