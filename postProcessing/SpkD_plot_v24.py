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
	fig=pylab.figure(1,figsize=(8,6))
	ax=fig.add_axes([0.13, 0.13, 0.7, 0.8],axisbg=bgColor)
	if LongSpikes and RepolSpikes and g['IncludeLongSpikes'].value:
		ax0=fig.add_axes([0.84, 0.13, 0.015, 0.8])
		ax1=fig.add_axes([0.855, 0.13, 0.015, 0.8])
	else:
		ax1=fig.add_axes([0.845, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
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
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.13, 0.13, 0.7, 0.8],axisbg=bgColor)
	if LongSpikes and RepolSpikes and g['IncludeLongSpikes'].value:
		ax0=fig.add_axes([0.84, 0.13, 0.015, 0.8])
		ax1=fig.add_axes([0.855, 0.13, 0.015, 0.8])
	else:
		ax1=fig.add_axes([0.845, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
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

def Densityplot(HdfFile, FileName, FieldName='', ProbFieldName=''\
, LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=1, cMax=10000, Res=8, logScale=True, cLabel=''\
, RepolSpikes=True, LongSpikes=True, ampThreshold=0, pThreshold=0, cThreshold=-1e12\
, nThreshold=0, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r, alpha=1.):
	assert cMin<cMax
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
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


def RateComparison(HdfFiles, FileName, RGBlist, FieldName=''\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=100, cMax=10000, Res=8, logScale=True, cLabel=''\
, RepolSpikes=True, LongSpikes=True, ampThreshold=0, pThreshold=0):
	assert cMin<cMax
	assert xMin<xMax
	X=np.zeros(((yMax-yMin)*Res,(xMax-xMin)*Res,3))
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8])
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
	for i in range(len(HdfFiles)):
		g=h5py.File(HdfFiles[i],'r')
		Sampling=g['Sampling'].value
		if 'CorrelationAnalysis/Probability' in g:
			p=g['CorrelationAnalysis/Probability'].value
		else:
			p=1.
		Amp=g['Amplitudes'].value
		pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
		Loc=g['Locations'].value[pAmpInd,:]
		cQb=False
		if not (FieldName==''):#plot the event property
			cQ=g[FieldName].value[pAmpInd]
			cQb=True
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
		if logScale:
			if cQb:
				X+=((np.log10(np.clip(Hq*1./np.clip(H,1,1e12),cMin,cMax))-np.log10(cMin))*1./(np.log10(cMax)-np.log10(cMin)))[:,:,None]*RGBlist[i,:][None,None,:]
			else:
				X+=((np.log10(np.clip(H*Res**2/4.,cMin,cMax))-np.log10(cMin))*1./(np.log10(cMax)-np.log10(cMin)))[:,:,None]*RGBlist[i,:][None,None,:]
		else:
			if cQb:
				X+=((np.clip(Hq*1./np.clip(H,1,1e12),cMin,cMax)-cMin)*1./(cMax-cMin))[:,:,None]*RGBlist[i,:][None,None,:]
			else:
				X+=((np.clip(H*Res**2/4.,cMin,cMax)-cMin)*1./(cMax-cMin))[:,:,None]*RGBlist[i,:][None,None,:]
	ax1.imshow(RGBlist[:,None,:]*np.ones(5)[None,:,None],origin='lower',aspect='auto',interpolation='none')
	ax1.yaxis.set_ticks_position('right')
	ax1.yaxis.set_label_position('right')
	ax1.set_xticks(np.array([]))
	ax1.set_yticks(np.arange(len(HdfFiles)))
	ax1.set_yticklabels(np.arange(len(HdfFiles),dtype=int))
	if not (cLabel==''):
		ax1.set_ylabel(cLabel,fontsize='small')
	else:
		ax1.set_ylabel('Phase',fontsize='small')
	ax.imshow(np.clip(X,0,1), aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return

def Matrixplot(HdfFile, FileName, FieldName='CorrelationAnalysis/Noise'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=1., Res=3, logScale=False, cLabel=''\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	assert cMin<cMax
	assert xMin<xMax
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
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
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
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
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.15, 0.13, 0.7, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.865, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
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

def Shapesplot_stats(HdfFile, FileName, FieldName='Shapes', ProbFieldName='CorrelationAnalysis/Probability'\
	, LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8\
, ampThreshold=0, pThreshold=0, cThreshold=-1e12\
, nThreshold=10, alpha=0.5):
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(8,4))
	ax=fig.add_axes([0.07, 0.13, 0.4, 0.8])
	ax0=fig.add_axes([0.54, 0.13, 0.4, 0.8])
	ax.tick_params(pad=4)
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	if not (ProbFieldName==''):
		p=g[ProbFieldName].value
	else:
		p=1.
	Amp=g[AmpFieldName].value
	pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	Loc=g[LocFieldName].value[pAmpInd,:]
	Shapes=g['Shapes'].value[pAmpInd,:]
	print Shapes.shape, Loc.shape
	sInd=g['RepolarizingSpikes'].value[pAmpInd]
	H=np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	H=H[1:,1:]+H[:-1,1:]+H[:-1,:-1]+H[1:,:-1]
	Hq1=np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res)\
	,weights=np.max(Shapes[sInd,:8],axis=1))[0]
	Hq2=np.clip(np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res)\
	,weights=np.clip(-np.min(Shapes[sInd,5:10],axis=1),4,1e12))[0],5,1e12)
	Hq3=np.histogram2d(Loc[sInd,0],Loc[sInd,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res)\
	,weights=np.max(Shapes[sInd,7:18],axis=1))[0]
	Hq1=Hq1[1:,1:]+Hq1[:-1,1:]+Hq1[:-1,:-1]+Hq1[1:,:-1]
	Hq2=Hq2[1:,1:]+Hq2[:-1,1:]+Hq2[:-1,:-1]+Hq2[1:,:-1]
	Hq3=Hq3[1:,1:]+Hq3[:-1,1:]+Hq3[:-1,:-1]+Hq3[1:,:-1]
	Hx=np.concatenate(((2*Hq1/Hq2)[:,:,None]\
	,(((Hq2*1./np.clip(H,nThreshold,1e12))-4)/4.)[:,:,None],(2*Hq3/Hq2)[:,:,None]),axis=2)
	ax.imshow(np.clip(Hx,0,1), aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax0.plot(0.5*Hx[:,:,0].flatten(),0.5*Hx[:,:,2].flatten(),'k,',alpha=alpha)
	ax0.set_ylim(0,0.5)
	ax0.set_xlim(0,0.5)
	ax0.set_xlabel('pre-peak/amplitude')
	ax0.set_ylabel('post-peak/amplitude')
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return

def GlobalVoltageplot(HdfFile, FileName, cMin=-0.2, cMax=0.3, CorrCmap=pylab.cm.hsv):
	fig=pylab.figure(1,figsize=(16,8))
	ax=fig.add_axes([0.1, 0.45, 0.35, 0.5])
	ax1=fig.add_axes([0.88, 0.4, 0.015, 0.55])
	#ax2=fig.add_axes([0.0625, 0.1, 0.125, 0.25])
	#ax2=fig.add_axes([0.5, 0.4, 0.125, 0.25])
	#ax2=fig.add_axes([0.5, 0.7, 0.125, 0.25])
	xPos=np.array([0.5,0.625,0.75,0.5,0.625,0.75,0.0625,0.1875,0.3125,0.4375,0.5625,0.6875,0.8125])
	yPos=np.array([0.7,0.7,0.7,0.4,0.4,0.4,0.1,0.1,0.1,0.1,0.1,0.1,0.1])
	ax.tick_params(pad=8)
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

def SurrogateStatsplot(HdfFile, FileName, tMin=0, tMax=20, Cmap=pylab.cm.jet):
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	SourceHist=g['SurrogateData/SourceHist'].value
	RateDependency=g['SurrogateData/RateDependency'].value
	FCluster=g['SurrogateData/FCluster'].value
	RefInd=g['SurrogateData/RefInd'].value
	OldInd=g['SurrogateData/OldInd'].value
	NewTimes=g['SurrogateData/NewTimes'].value
	NewIndMap=g['SurrogateData/NewIndMap'].value
	NewLoc=g['Locations'].value[NewIndMap]
	pylab.figure(1,figsize=(16,12))
	pylab.subplot(221)
	pylab.semilogy(np.cumsum(SourceHist[::-1])[::-1][:60],'b-')
	pylab.semilogy(SourceHist[:60],'k-')
	pylab.ylabel('count (blue: cumulative)')
	pylab.xlabel('copy frequency')
	pylab.subplot(222)
	pylab.imshow((np.cumsum(RateDependency[:,:60],axis=0)*1./np.clip(SourceHist,1,1e12)[None,:60])\
	,cmap=Cmap,aspect='auto',interpolation='none',origin='lower',extent=(0,60,0,1))
	pylab.xlabel('copy frequency')
	pylab.ylabel('time/100 events/sec')
	cbar=pylab.colorbar()
	cbar.set_label('cumulative density')
	pylab.subplot(223)
	Colors=np.argsort(scipy.rand(FCluster.max()+1))
	CM=np.zeros(192**2)-1
	CM[RefInd]=Colors[FCluster]
	pylab.imshow(np.reshape(CM,(192,192)),cmap=Cmap,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,64,0,64))
	cbar=pylab.colorbar()
	cbar.set_label('cluster id')
	pylab.subplot(224)
	Ind=(NewTimes>tMin*Sampling)*(NewTimes<tMax*Sampling)
	pylab.plot(NewTimes[Ind]*1./Sampling,NewLoc[Ind,0], 'r,', alpha=0.1)
	Ind=(NewTimes[OldInd]>tMin*Sampling)*(NewTimes[OldInd]<tMax*Sampling)
	pylab.plot(NewTimes[OldInd[Ind]]*1./Sampling,NewLoc[OldInd[Ind],0], 'k,', alpha=0.2)
	pylab.xlabel('time')
	pylab.ylim(0,64)
	#pylab.ylabel('Location(y)')
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return


def PrintStats(HdfFile,StatsFile, onlineVersion=False):
	if onlineVersion:
		HighThr=8
	else:
		HighThr=2.
	g=h5py.File(HdfFile,'r')
	f=file(StatsFile,'w')
	f.write('Sampling rate: %i \n'%(g['Sampling'].value))
	f.write('Duration of recording: %i s \n'%(g['tMax'].value))
	f.write('Number of recorded channels: %i \n \n'%(len(g['recordedChannels'].value)))
	if not onlineVersion:
		f.write('Amplitude threshold: %i\n'%(g['AmplitudeThresholds'].value[0]))
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
	f.write('high amplitude (%2.1f) stats: \n'%HighThr)
	f.write('Number of true spikes %i \n'%(np.sum((g['CorrelationAnalysis/Probability'].value\
	[(g['Amplitudes'].value>=HighThr)]))))
	f.write('Fraction of noisy spikes %4.3f \n \n'%(1.-np.mean((g['CorrelationAnalysis/Probability'].value\
	[(g['Amplitudes'].value>=HighThr)]))))
	if not onlineVersion:
		f.write('Number of noisy spikes (repolarizing) %i \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,0])))
		f.write('Number of true spikes (repolarizing) %i \n'%(np.sum(((1.-g['CorrelationAnalysis/Noise'].value)*NSpikes)[:,0])))
		f.write('Fraction of noisy spikes (repolarizing) %4.3f \n \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,0])\
		*1./np.sum(NSpikes[:,0])))
		f.write('Number of noisy long spikes %i \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,1])))
		f.write('Fraction of noisy long spikes %4.3f \n \n'%(np.sum((g['CorrelationAnalysis/Noise'].value*NSpikes)[:,1])\
		*1./np.sum(NSpikes[:,1])))
		f.write('high amplitude (%2.1f) stats: \n'%HighThr)
		f.write('Number of true spikes (repolarizing) %i \n'%(np.sum((g['CorrelationAnalysis/Probability'].value\
		[(g['Amplitudes'].value>=HighThr)*g['RepolarizingSpikes'].value]))))
		f.write('Fraction of noisy spikes (repolarizing) %4.3f \n \n'%(1.-np.mean((g['CorrelationAnalysis/Probability'].value\
		[(g['Amplitudes'].value>=HighThr)*g['RepolarizingSpikes'].value]))))
		f.write('Number of true long spikes %i \n'%(np.sum((g['CorrelationAnalysis/Probability'].value\
		[(g['Amplitudes'].value>=HighThr)*(True-g['RepolarizingSpikes'].value)]))))
		f.write('Fraction of noisy long spikes %4.3f \n \n'%(1.-np.mean((g['CorrelationAnalysis/Probability'].value\
		[(g['Amplitudes'].value>=HighThr)*(True-g['RepolarizingSpikes'].value)]))))
		f.write('Cluster \n')
		f.write('Number of clusters %i \n'%(len(g['Cluster/NCount'].value)-1))
		f.write('Number of events in clusters %i \n'%(np.sum(g['Cluster/NCount'].value[1:])))
	f.close()
	g.close()
	return

def Excitabilityplot(HdfFile, FileName, FieldName='Excitability3/relE', ProbFieldName=''\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, cMin=0, cMax=3, Res=8, cLabel=''\
, ampThreshold=0, pThreshold=0, cThreshold=-1e12\
, nThreshold=0, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r, alpha=1.):
	assert cMin<cMax
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.41, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.52, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile,'r')
	Sampling=g['Sampling'].value
	useInd=False
	if not (IndFieldName==''):
		Ind=g[IndFieldName].value
		useInd=True
	if not (ProbFieldName==''):
		p=g[ProbFieldName].value
		if useInd:
			p=p[Ind]
	else:
		p=1.
	Amp=g[AmpFieldName].value
	if useInd:
		Amp=Amp[Ind]
	if not (FieldName==''):#plot the event property
		cQ=g[FieldName].value
		if useInd:
			cQ=cQ[Ind]
		pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)*(cQ>=cThreshold)
		cQ=cQ[pAmpInd]
	else:
		pAmpInd=(p>=pThreshold)*(Amp>=ampThreshold)
	if useInd:
		Loc=g[LocFieldName].value[Ind,:][pAmpInd,:]
	else:
		Loc=g[LocFieldName].value[pAmpInd,:]
	Amp=Amp[pAmpInd]
	A=100.-(np.argsort(np.argsort(Amp))*100./len(Amp))
	L=len(A)
	#Havg=np.zeros(((yMax-yMin)*Res,(yMax-yMin)*Res,100))
	#Hvar=np.zeros(((yMax-yMin)*Res,(yMax-yMin)*Res,100))
	#Hn=np.zeros(((yMax-yMin)*Res,(yMax-yMin)*Res,100))
	H=np.histogramdd(np.concatenate((Loc,A[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H=H[1:,1:,:]+H[:-1,1:,:]+H[:-1,:-1,:]+H[1:,:-1,:]
	Hq=np.histogramdd(np.concatenate((Loc,A[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)),weights=cQ)[0]
	Hs=np.histogramdd(np.concatenate((Loc,A[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)),weights=cQ**2)[0]
	Hq=Hq[1:,1:,:]+Hq[:-1,1:,:]+Hq[:-1,:-1,:]+Hq[1:,:-1,:]
	Hs=Hs[1:,1:,:]+Hs[:-1,1:,:]+Hs[:-1,:-1,:]+Hs[1:,:-1,:]
	H=np.clip(np.cumsum(H,axis=2),2,1e12)
	Hq=np.cumsum(Hq,axis=2)*1./H
	Hs=(np.cumsum(Hs,axis=2)*1./H-Hq**2)*H*1./(H-1)+1e12*(H<4)
	fHighAmp=99-np.argmin((Hs*1./H)[:,:,::-1],axis=2)
	IndX=fHighAmp.flatten()
	print IndX.shape, IndX.min(), IndX.mean(), np.percentile(IndX,80)
	Havg=np.reshape(np.reshape(Hq,(-1,100))[np.arange(IndX.shape[0]),IndX],((yMax-yMin)*Res,-1))
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	if not (cLabel==''):
		cbar1.set_label(cLabel,fontsize='small')
	elif not (FieldName==''):
		cbar1.set_label(FieldName,fontsize='small')
	else:
		cbar1.set_label('Event count',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=SpikesCmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('fraction of events',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	kMax=int(np.ceil(1./alpha))
	if (kMax>1):
		for k in range(kMax):
			ax.imshow(np.ma.array(Havg,mask=True-((H[:,:,-1]>=(2**(kMax-k-1)))*(H[:,:,-1]>nThreshold)))\
			,cmap=SpikesCmap, aspect='equal'\
			,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax), alpha=alpha)
			ax0=fig.add_axes([0.845, 0.13+(0.025*k)/kMax, 0.025/kMax, 0.8])
			ax0.set_frame_on(False)
			ax0.set_yticks(np.array([]))
			ax0.set_xticks(np.array([]))
			ax0.imshow(np.outer(np.ones(100),np.ones(5)),cmap=pylab.cm.gray,origin='lower',aspect='auto')
			ax0.imshow(np.outer(np.linspace(0,1,100),np.ones(5)),alpha=alpha\
			,cmap=LongSpikesCmap,origin='lower',aspect='auto')
	else:#
		ax.imshow(np.ma.array(Hq[:,:,-1]-Havg,mask=True-((H[:,:,-1]>=1)*(H[:,:,-1]>nThreshold)))\
		,vmin=cMin, vmax=cMax,cmap=SpikesCmap, aspect='equal'\
		,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax), alpha=alpha)
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(np.ma.array(fHighAmp/100.,mask=True-((H[:,:,-1]>=1)*(H[:,:,-1]>nThreshold)))\
	,vmin=0., vmax=1., cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return

def IVarPlot(HdfFile, FileName):
	#f=h5py.File(NoisyChFile,'r+')
	#NoisyAreas=f['NoisyAreas'].value
	#f.close()
	#nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	#Amp=g['Amplitudes'].value
	Times=(g['Times'].value)[SingleSpk]
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res=6
	nBins=64*Res
	Ns=2
	tau=1./Sampling
	Spikes0=np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)+Ns
	Spikes1=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)+Ns
	N=Spikes0.shape[0]
	Last=np.zeros((nBins+2*Ns,nBins+2*Ns))
	dt=[]
	for i in range(N):
		j=Spikes0[i]
		k=Spikes1[i]
		dt.append(Times[i]-np.max(Last[j-Ns:j+Ns+1,k-Ns:k+Ns+1]))
		Last[j,k]=Times[i]
	dtE=np.exp(-np.array(dt)*tau)
	dtN=np.histogram2d(Spikes0,Spikes1,bins=(np.arange(Ns,nBins+Ns),np.arange(Ns,nBins+Ns)))[0]
	dtAvg=np.histogram2d(Spikes0,Spikes1,bins=(np.arange(Ns,nBins+Ns),np.arange(Ns,nBins+Ns)),weights=dtE)[0]
	dtA=dtAvg*1./np.clip(dtN,1,1e12)
	dtSq=np.sqrt(np.histogram2d(Spikes0,Spikes1,bins=(np.arange(Ns,nBins+Ns),np.arange(Ns,nBins+Ns))\
	,weights=((dtE-dtA[Spikes0-Ns,Spikes1-Ns])**2))[0])
	Istd=(dtSq/np.clip(dtN,1,1e12))*1./np.clip(1.-dtA,1e-6,1e12)
	print dtA.min(),dtA.max(),dtA.mean()
	print Istd.min(),Istd.max(),Istd.mean()
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8])
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=pylab.cm.spectral,norm=mplt.colors.Normalize(vmin=0,vmax=0.5))
	cbar1.set_label('ISI variability',fontsize='small')
	ax.imshow(np.ma.array(Istd,mask=(dtN<5)),vmin=0., vmax=0.5\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,64,0,64))
	ax.set_ylim(0,64)
	ax.set_xlim(0,64)
	ax.set_xticks(np.array([0,32,64]))
	ax.set_yticks(np.array([0,32,64]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return


def IslowVarPlot(HdfFile, FileName):
	#f=h5py.File(NoisyChFile,'r+')
	#NoisyAreas=f['NoisyAreas'].value
	#f.close()
	#nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=g['Amplitudes'].value[SingleSpk]
	Times=(g['Times'].value)[SingleSpk]
	Loc=Loc[Amp>4.,:]
	Times=Times[Amp>4.]
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res=6
	nBins=64*Res
	Ns=2
	tau=0.01/Sampling
	Spikes0=np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)+Ns
	Spikes1=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)+Ns
	N=Spikes0.shape[0]
	Last=np.zeros((nBins+2*Ns,nBins+2*Ns))
	dt=[]
	for i in range(N):
		j=Spikes0[i]
		k=Spikes1[i]
		dt.append(Times[i]-np.max(Last[j-Ns:j+Ns+1,k-Ns:k+Ns+1]))
		Last[j,k]=Times[i]+ 1./tau*np.clip(np.exp(-dt[-1]*tau),0,dt[-1]*tau)
	dtN=np.histogram2d(Spikes0,Spikes1,bins=(np.arange(0,nBins+2*Ns+1),np.arange(0,nBins+2*Ns+1)))[0]
	dtN0=np.zeros(dtN.shape)
	for j in range(2*Ns+1):
		for k in range(2*Ns+1):
			dtN0[Ns:-Ns,Ns:-Ns]+=dtN[j:nBins+j,k:nBins+k]
	dtN=dtN[Ns:-Ns,Ns:-Ns]
	dtN0=dtN0[Ns:-Ns,Ns:-Ns]
	dtE=np.exp(-np.array(dt)*tau)
	dtAvg=np.histogram2d(Spikes0,Spikes1,bins=(np.arange(Ns,nBins+Ns+1),np.arange(Ns,nBins+Ns+1)),weights=dtE)[0]
	dtA=dtAvg*1./np.clip(dtN,1,1e12)
	dtSq=np.sqrt(np.histogram2d(Spikes0,Spikes1,bins=(np.arange(Ns,nBins+Ns+1),np.arange(Ns,nBins+Ns+1))\
	,weights=((dtE-dtA[Spikes0-Ns,Spikes1-Ns])**2))[0])
	Istd=(dtSq/np.clip(dtN,1,1e12)*np.clip(np.exp(dtN0/np.mean(dtN0[dtN0>=1000])-1.),0,1))#/np.sqrt(np.clip(dtA,1e-6,1e12))
	print dtA.min(),dtA.max(),dtA.mean()
	print Istd.min(),Istd.max(),Istd.mean()
	fig=pylab.figure(1,figsize=(8,7))
	ax=fig.add_axes([0.07, 0.13, 0.7, 0.8])
	ax1=fig.add_axes([0.785, 0.13, 0.025, 0.8])
	ax.tick_params(pad=8)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=pylab.cm.spectral,norm=mplt.colors.Normalize(vmin=0,vmax=0.1))
	cbar1.set_label('rate variability',fontsize='small')
	ax.imshow(np.ma.array(Istd,mask=(dtN<5)),vmin=0., vmax=0.1\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,64,0,64))
	ax.set_ylim(0,64)
	ax.set_xlim(0,64)
	ax.set_xticks(np.array([0,32,64]))
	ax.set_yticks(np.array([0,32,64]))
	pylab.savefig(FileName)
	pylab.close()
	g.close()
	return
	
def FFPlot(HdfFile, FileName):
	#f=h5py.File(NoisyChFile,'r+')
	#NoisyAreas=f['NoisyAreas'].value
	#f.close()
	#nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=g['Amplitudes'].value[SingleSpk]
	Times=np.array((g['Times'].value)[SingleSpk],dtype=int)
	Loc=Loc[Amp>3.,:]
	Times=Times[Amp>3.]
	Times=np.concatenate((Times,np.array([tMax*Sampling+1])))
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res=5
	Ns0=2
	Ns=np.array([2])
	nBins=64*Res+2*Ns
	tau=0.01/Sampling
	Spikes0=np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)+Ns
	Spikes1=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)+Ns
	N=Spikes0.shape[0]
	#Spikes=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)+Ns0\
	#+(np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)+Ns0)*(nBins)
	#N=Spikes.shape[0]
	i=0
	for r in range(1):
		Last=[]
		FF=[]
		FFN=[]
		for i in range(14):
			Last.append(np.zeros((nBins+2*Ns,nBins+2*Ns,2),dtype=int))
			FF.append(np.zeros((nBins,nBins),dtype=int))
			FFN.append(np.zeros((nBins,nBins),dtype=int))
			Nff=np.zeros(14)
		for t in range(int(tMax*Sampling)/10):
			Last[0][:,:,t%2]=0
			if t%10000==0:
				print t*10./Sampling
			while Times[i]<=10*t:
				j=Spikes0[i]
				k=Spikes1[i]
				Last[0][j-Ns:j+Ns+1,k-Ns:k+Ns+1,t%2]+=1
				i+=1
			for k in range(1,14):
				if ((t+1)%(2**(k-1))==0) and t>=(2**k):
					Last[k][:,:,(t+1)/(2**(k-1))%2]=np.sum(Last[k-1],axis=2)
					FF[k]+=(Last[k][Ns:-Ns,Ns:-Ns,(t+1)/(2**(k-1))%2])**2
					FFN[k]+=(Last[k][Ns:-Ns,Ns:-Ns,(t+1)/(2**(k-1))%2])
					Nff[k]+=1
			FF[0]+=(Last[0][Ns:-Ns,Ns:-Ns,(t+1)%2])**2
			FFN[0]+=(Last[0][Ns:-Ns,Ns:-Ns,(t+1)%2])
			Nff[0]+=1
	#rate histogram
	dtN=np.histogram2d(Spikes0,Spikes1,bins=(np.arange(0,nBins+2*Ns+1),np.arange(0,nBins+2*Ns+1)))[0]
	dtN0=np.zeros(dtN.shape)
	for j in range(2*Ns+1):
		for k in range(2*Ns+1):
			dtN0[Ns:-Ns,Ns:-Ns]+=dtN[j:nBins+j,k:nBins+k]
	dtN=dtN[Ns:-Ns,Ns:-Ns]
	dtN0=dtN0[Ns:-Ns,Ns:-Ns]
	#expected count
	#Eff=dtN0[:,:,None]*(2.**(np.arange(16))/(tMax*Sampling))[None,None,:]
	#normalized FF
	FF0=np.zeros((nBins,nBins,14))
	for i in range(14):
		FF0[:,:,i]=(FF[i]-FFN[i]**2*1./np.clip(Nff[i],1,1e12))*1./np.clip(FFN[i],1,1e12)
	#what to visualize?
	# -- variance for each timescale (densityplot)
	D=np.zeros((90,14))
	Savg=np.zeros(14)
	Tavg=np.zeros((nBins,nBins))
	Tasy=np.zeros((nBins,nBins))
	for i in range(14):
		D[:,i]=np.histogram(FF0[:,:,i],bins=np.logspace(-1,2,91))[0]
		Savg[i]=np.exp(np.mean(np.log(FF0[:,:,i][dtN0>100])))
		Tavg[dtN0>10]+=np.log10(FF0[:,:,i][dtN0>10]/Savg[i])
		Tasy[dtN0>10]+=((i-7.5)/32.)*np.log10(FF0[:,:,i][dtN0>10]/Savg[i])
	Tavg=Tavg/14.
	# -- difference to local average?
	# -- average
	# -- temporal asymmetry (first 8 - last 8, or better center of mass?)
	fig3=pylab.figure(3,figsize=(6.4,4.8))
	#pylab.subplots_adjust(left=0.12, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3)
	ax10=fig3.add_axes([0.05, 0.53, 0.3, 0.38])
	ax11=fig3.add_axes([0.55, 0.53, 0.3, 0.38])
	ax12=fig3.add_axes([0.05, 0.04, 0.3, 0.38])
	ax13=fig3.add_axes([0.55, 0.04, 0.3, 0.38])
	ax14=fig3.add_axes([0.4, 0.53, 0.015, 0.38])
	ax15=fig3.add_axes([0.4, 0.04, 0.015, 0.38])
	ax16=fig3.add_axes([0.9, 0.53, 0.015, 0.38])
	ax17=fig3.add_axes([0.9, 0.04, 0.015, 0.38])
	ax10.tick_params(pad=4)
	ax11.tick_params(pad=4)
	ax12.tick_params(pad=4)
	ax13.tick_params(pad=4)
	fig3.text(0.01,0.93,'A',fontsize='xx-large')
	fig3.text(0.51,0.93,'B',fontsize='xx-large')
	fig3.text(0.01,0.45,'C',fontsize='xx-large')
	fig3.text(0.51,0.45,'D',fontsize='xx-large')
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax14,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(1,1e5))
	cbar0.set_label(r'count',fontsize='x-small')
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax15,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(0.1,1e1))
	cbar0.set_label(r'avg. mult. deviation',fontsize='x-small')
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax16,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(0.1,1e1))
	cbar0.set_label(r'rel.slope',fontsize='x-small')
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax17,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(0.1,1e2))
	cbar0.set_label(r'Fano Factor',fontsize='x-small')
	
	ax10.imshow(np.log10(np.clip(D,1,1e12)),vmin=0., vmax=5.\
	,cmap=pylab.cm.spectral, aspect='auto',interpolation='none',origin='lower',extent=(-0.5,13.5,0,30))
	ax11.imshow(np.ma.array(Tavg,mask=(dtN0<10)),vmin=-1., vmax=1.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,64,0,64))
	ax12.imshow(np.ma.array(Tasy,mask=(dtN0<10)),vmin=-1., vmax=1.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,64,0,64))
	ax13.imshow(np.ma.array(np.log10(FF0[:,:,13]),mask=(dtN0<10)),vmin=-1., vmax=2.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,64,0,64))
	ax10.set_ylim(0,30)
	ax10.set_xlim(-0.5,13.5)
	ax10.set_yticks(np.array([0,10,20,30]))
	ax10.set_yticklabels(np.array([0.1,1,10,100]))
	ax10.set_xticks(np.array([0,6,13]))
	ax10.set_xticklabels(np.array(['%1.3f'%(10./Sampling),'%1.2f'%(640./Sampling),'%2.1f'%(81920./Sampling)]))
	ax10.set_xlabel('T/s')
	ax10.set_xlabel('Fano Factor')
	ax11.set_ylim(0,64)
	ax11.set_xlim(0,64)
	ax11.set_xticks(np.array([0,32,64]))
	ax11.set_yticks(np.array([0,32,64]))
	ax12.set_ylim(0,64)
	ax12.set_xlim(0,64)
	ax12.set_xticks(np.array([0,32,64]))
	ax12.set_yticks(np.array([0,32,64]))
	ax13.set_ylim(0,64)
	ax13.set_xlim(0,64)
	ax13.set_xticks(np.array([0,32,64]))
	ax13.set_yticks(np.array([0,32,64]))
	pylab.savefig(FileName)
	np.save(FileName,FF0)
	pylab.close()
	g.close()
	return


def FFxtPlot(HdfFile, FileName):
	#f=h5py.File(NoisyChFile,'r+')
	#NoisyAreas=f['NoisyAreas'].value
	#f.close()
	#nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=g['Amplitudes'].value[SingleSpk]
	Times=np.array((g['Times'].value)[SingleSpk],dtype=int)
	#Loc=Loc[Amp>3.,:]
	#Times=Times[Amp>3.]
	#Times=np.concatenate((Times,np.array([tMax*Sampling+1])))
	I=np.concatenate((np.array([0]),np.cumsum(np.histogram(Times,bins=np.arange(0,int(tMax*Sampling)+20,20))[0])))
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res=8
	#Ns=3###use periodic boundaries
	#Ns=np.array([2])
	nBins=64*Res
	tau=0.01/Sampling
	Spikes=np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)*nBins\
	+np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)
	#print np.max(Spikes%nBins)
	#N=Spikes0.shape[0]
	#Spikes=np.clip(np.array(Res*(Loc[:,1]),dtype=int),0,nBins-1)+Ns0\
	#+(np.clip(np.array(Res*(Loc[:,0]),dtype=int),0,nBins-1)+Ns0)*(nBins)
	#N=Spikes.shape[0]
	#i=0
	for r in range(1):
		Last=[]
		FF=[]
		FFN=[]
		#nB=nBins-3
		for i in range(14):
			nB=nBins-3
			Last.append(np.zeros((nB**2,2),dtype=int))
			FF.append(np.zeros((nB)**2,dtype=int))#half electrode...
			FFN.append(np.zeros((nB)**2,dtype=int))
			for j in range(5):
				nB=(nBins/2**(j+1)-3)
				Last.append(np.zeros((nB**2,2),dtype=int))
				FF.append(np.zeros((nB)**2,dtype=int))
				FFN.append(np.zeros((nB)**2,dtype=int))
			Nff=np.zeros(14)
		for t in range(int(tMax*Sampling)/20):
			if I[t]<I[t+1]:
				#print Spikes[I[t]:I[t+1]]*1./nBins
				X=np.reshape(np.bincount(Spikes[I[t]:I[t+1]],minlength=nBins**2),(nBins,nBins))
				#print np.nonzero(X)
			else:
				X=np.zeros((nBins,nBins),dtype=int)
			Y=X.copy()#X[3:,3:]#(nBins-3)
			for ky in range(1,4):
				Y[3:,3:]+=X[3-ky:-ky,3:]
			for kx in range(1,4):
				Y[3:,3:]+=X[3:,3-kx:-kx]
				for ky in range(1,4):
					Y[3:,3:]+=X[3-ky:-ky,3-kx:-kx]
			Last[0][:,t%2]=Y[3:,3:].flatten()
			Y=np.reshape(Y,(4,Y.shape[0]/4,4,Y.shape[1]/4))
			for j in range(1,6):
				X=Y.copy()
				Y=X[:,::2,:,::2]
				Y+=X[:,1::2,:,::2]
				Y+=X[:,::2,:,1::2]
				Y+=X[:,1::2,:,1::2]
				Last[j][:,t%2]=(np.reshape(Y,(Y.shape[1]*4,Y.shape[3]*4))[3:,3:]).flatten()
				#print np.nonzero(Last[j][:,t%2])
			if t%10000==0:
				print t*20./Sampling
			for k in range(1,14):
				if ((t+1)%(2**(k-1))==0) and t>=(2**k):
					Nff[k]+=1
					for j in range(6):
						#print k,j
						Last[6*k+j][:,(t+1)/(2**(k-1))%2]=np.sum(Last[6*(k-1)+j],axis=1)
						FF[6*k+j]+=(Last[6*k+j][:,(t+1)/(2**(k-1))%2])**2
						FFN[6*k+j]+=(Last[6*k+j][:,(t+1)/(2**(k-1))%2])
			for j in range(6):
				FF[j]+=(Last[j][:,(t+1)%2])**2
				FFN[j]+=(Last[j][:,(t+1)%2])
			Nff[0]+=1
	#normalized FF
	#FF0=np.zeros((nBins,nBins,14))
	#for i in range(14):
	#	FF0[:,:,i]=(FF[i]-FFN[i]**2*1./np.clip(Nff[i]-1,1,1e12))*1./np.clip(FFN[i],1,1e12)
	#save to .hdf
	f=h5py.File(FileName+'.hdf','w')
	f.create_dataset('Nff', data=Nff)
	for i in range(14):
		for j in range(6):
			f.create_dataset('Sq%i_%i'%(i,j), data=FF[i*6+j])
			f.create_dataset('N%i_%i'%(i,j), data=FFN[i*6+j])
	f.close()
	g.close()
	return

def FFxtPlot2(HdfFile, FileName):
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	g.close()
	g=h5py.File(FileName+'.hdf','r+')
	Nff=g['Nff'].value
	Sq=[]
	N=[]
	for i in range(14):
		for j in range(6):
			Sq.append(g['Sq%i_%i'%(i,j)].value)
			N.append(g['N%i_%i'%(i,j)].value)
	g.close()
	Res=8
	nBins=64*Res
	Area=4.**np.arange(-1,5)# in number of electrodes (j)
	T=20./Sampling*2.**np.arange(14)
	Tlabel=np.array(['%1.4f'%T[1],'%1.3f'%T[4],'%1.2f'%T[8],'%2.1f'%T[12]])
	Alabel=np.array(['0.25','4','64'])
	NbyA=[]
	FF0=[]
	nB=[]
	TbyA=np.zeros(6*14)
	for i in range(14):
		#nB.append(nBins-3)
		for j in range(6):
			#if j>0:
			nB.append(nBins/2**j-3)
			#NbyA.append((3*np.argsort(np.argsort(N[i*6+j])))/np.size(N[i*6+j]))#classify into 3 groups
			TbyA[i*6+j]=(i-j+5)#colourscale[0:24)#T[i]/A[j]#has discrete values...use log and colourscale?
			FF0.append((Sq[i*6+j]-N[i*6+j]**2*1./np.clip(Nff[i],1,1e12))*1./np.clip(N[i*6+j],1,1e12))
			#NbyA.append((3*np.argsort(np.argsort(FF0[i*6+j])))/np.size(FF0[i*6+j]))#classify into 3 groups
	fig0=pylab.figure(3,figsize=(9,6))
	#pylab.subplots_adjust(left=0.12, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3)
	I=np.array([1,4,8,12])
	J=np.array([0,2,4])
	print N
	#for i in range(30):
	#	print FF0[i].min(), FF0[i].max(), NbyA[i].mean(), Sq[i].mean(), N[i].mean()
	for i in range(4):
		for j in range(3):
			ax=fig0.add_axes([0.05+0.2*i,0.05+0.3*j,0.2,0.3])
			ax.tick_params(pad=4)
			ax.imshow(np.log10(np.reshape(np.ma.array(FF0[I[i]*6+J[j]],mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]]))),vmin=-1., vmax=3.\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			ax.set_xticks(np.array([]))
			ax.set_yticks(np.array([]))
			if i==0:
				ax.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax.set_xlabel(Tlabel[i]+' s')
	ax0=fig0.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(0.1,1e3))
	cbar0.set_label(r'Fano Factor',fontsize='x-small')
	pylab.savefig(FileName+'_spatial.png')
	pylab.close()
	
	fig0=pylab.figure(3,figsize=(8,5))
	ax0=fig0.add_axes([0.1,0.58,0.25,0.4])
	ax1=fig0.add_axes([0.37,0.58,0.25,0.4])
	ax2=fig0.add_axes([0.64,0.58,0.25,0.4])
	ax3=fig0.add_axes([0.9,0.58,0.025,0.4])
	ax4=fig0.add_axes([0.1,0.09,0.25,0.4])
	ax5=fig0.add_axes([0.37,0.09,0.25,0.4])
	ax6=fig0.add_axes([0.64,0.09,0.25,0.4])
	ax7=fig0.add_axes([0.9,0.09,0.025,0.4])
	
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=pylab.cm.spectral\
	,norm=mplt.colors.LogNorm(0.5,500))
	cbar0.set_label(r'Fano Factor',fontsize='x-small')
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=pylab.cm.spectral\
	,norm=mplt.colors.LogNorm(20./Sampling/2.**4,20./Sampling*2.**14))
	cbar0.set_label(r'(T/s)/sqrt(area/electrodes)',fontsize='x-small')
	FFhist=np.zeros((19,20,3))
	FFhistN=np.zeros((19,20,3))
	for i in range(14):
		for j in range(6):
			#for k in range(3):
			#	#if k==2:
			#	#	print np.min(FF0[i*6+j][(N[i*6+j]>=10)*(NbyA[i*6+j]==k)]),Nff[i],Sq[i*6+j].mean(), N[i*6+j].mean()
			FFhist[TbyA[i*6+j],:,j/2]+=np.histogram(N[i*6+j][(N[i*6+j]>=10)]*1./Nff[i]\
			,bins=(0.01*2.**np.arange(21)),weights=np.log10(FF0[i*6+j][(N[i*6+j]>=10)]))[0]
			FFhistN[TbyA[i*6+j],:,j/2]+=np.histogram(N[i*6+j][(N[i*6+j]>=10)]*1./Nff[i]\
			,bins=(0.01*2.**np.arange(21)))[0]
			Qa=np.argsort(N[i*6+j])
			Q=np.exp(np.diff(np.cumsum(np.concatenate((np.array([0]),\
			np.log(np.clip(N[i*6+j],10,1e12)/Nff[i]))))[::nB[i*6+j]])*1./nB[i*6+j])
			Qf=np.exp(np.diff(np.cumsum(np.concatenate((np.array([0]),\
			np.log(np.clip(FF0[i*6+j],0.01,1e4)))))[::nB[i*6+j]])*1./nB[i*6+j])
			if j<2:
				ax4.loglog(Q,Qf,'.'\
				,mew=0.001,mfc=pylab.cm.spectral(TbyA[i*6+j]/18.),ms=j/2.+1.)
			elif j<4:
				ax5.loglog(Q,Qf,'.'\
				,mew=0.001,mfc=pylab.cm.spectral(TbyA[i*6+j]/18.),ms=j/2.+1.)
			else:
				ax6.loglog(Q,Qf,'.'\
				,mew=0.001,mfc=pylab.cm.spectral(TbyA[i*6+j]/18.),ms=j/2.+1.)
	ax4.set_ylabel('Fano Factor')
	ax4.set_xlabel('log # events')
	ax4.set_ylim(0.1,100)
	ax5.set_ylim(1.,1000)
	ax6.set_ylim(10.,10000)
	ax4.set_xlim(1e-5,1e1)
	ax5.set_xlim(1e-4,1e2)
	ax6.set_xlim(1e-3,1e4)
	ax0.imshow(np.ma.array((FFhist[:,1:,0]+FFhist[:,:-1,0])*1./np.clip((FFhistN[:,1:,0]+FFhistN[:,:-1,0]),1,1e12)\
	,mask=((FFhistN[:,1:,0]+FFhistN[:,:-1,0])==0)),vmin=np.log10(0.5), vmax=np.log10(500.)\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,1,0,1))
	ax0.set_ylabel('log2((T/20 frames)-2log4(area/ch))')
	ax0.set_xlabel('log # events')
	ax1.imshow(np.ma.array((FFhist[:,1:,1]+FFhist[:,:-1,1])*1./np.clip((FFhistN[:,1:,1]+FFhistN[:,:-1,1]),1,1e12)\
	,mask=((FFhistN[:,1:,1]+FFhistN[:,:-1,1])==0)),vmin=np.log10(0.5), vmax=np.log10(500.)\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,1,0,1))
	ax2.imshow(np.ma.array((FFhist[:,1:,2]+FFhist[:,:-1,2])*1./np.clip((FFhistN[:,1:,2]+FFhistN[:,:-1,2]),1,1e12)\
	,mask=((FFhistN[:,1:,2]+FFhistN[:,:-1,2])==0)),vmin=np.log10(0.5), vmax=np.log10(500.)\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none',origin='lower',extent=(0,1,0,1))
	ax0.set_xticks(np.array([0,np.log2(1000.)/20.,np.log2(1000000.)/20.]))
	ax0.set_xticklabels(np.array(['1e-2','10','1e4']))
	ax1.set_xticks(np.array([0,np.log2(1000.)/20.,np.log2(1000000.)/20.]))
	ax1.set_xticklabels(np.array(['1e-2','10','1e4']))
	ax2.set_xticks(np.array([0,np.log2(1000.)/20.,np.log2(1000000.)/20.]))
	ax2.set_xticklabels(np.array(['1e-2','10','1e4']))
	ax0.set_yticks(np.array([1./36.,9./36,35./36.]))
	ax0.set_yticklabels(np.array(['-4','0','13']))
	ax1.set_yticks(np.array([]))
	ax1.set_yticklabels(np.array([]))
	ax2.set_yticks(np.array([]))
	ax2.set_yticklabels(np.array([]))
	pylab.savefig(FileName+'_scatter.png')
	pylab.close()
	return

def FFxtPlot3(HdfFile, FileName):
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	g.close()
	g=h5py.File(FileName+'.hdf','r+')
	Nff=g['Nff'].value
	Sq=[]
	N=[]
	for i in range(14):
		for j in range(6):
			Sq.append(g['Sq%i_%i'%(i,j)].value)
			N.append(g['N%i_%i'%(i,j)].value)
	g.close()
	Res=8
	nBins=64*Res
	Area=4.**np.arange(-1,5)# in number of electrodes (j)
	T=20./Sampling*2.**np.arange(14)
	Tlabel=np.array(['%1.4f'%T[8],'%1.3f'%T[9],'%1.2f'%T[10],'%2.1f'%T[11]])
	Alabel=np.array(['0.25','1','4'])
	NbyA=[]
	FF0=[]
	nB=[]
	TbyA=np.zeros(6*14)
	for i in range(14):
		#nB.append(nBins-3)
		for j in range(6):
			#if j>0:
			nB.append(nBins/2**j-3)
			#NbyA.append((3*np.argsort(np.argsort(N[i*6+j])))/np.size(N[i*6+j]))#classify into 3 groups
			TbyA[i*6+j]=(i-j+5)#colourscale[0:24)#T[i]/A[j]#has discrete values...use log and colourscale?
			FF0.append(np.clip((Sq[i*6+j]-N[i*6+j]**2*1./np.clip(Nff[i],1,1e12))*1./np.clip(N[i*6+j],1,1e12),1e-6,1e12))
			#NbyA.append((3*np.argsort(np.argsort(FF0[i*6+j])))/np.size(FF0[i*6+j]))#classify into 3 groups
	fig0=pylab.figure(1,figsize=(9,6))
	fig2=pylab.figure(2,figsize=(9,6))
	fig3=pylab.figure(3,figsize=(9,6))
	#pylab.subplots_adjust(left=0.12, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3)
	I=np.array([8,9,10,11])
	J=np.array([0,1,2])
	J2=np.array([1,2,3])
	#print N
	#for i in range(30):
	#	print FF0[i].min(), FF0[i].max(), NbyA[i].mean(), Sq[i].mean(), N[i].mean()
	for i in range(4):
		for j in range(3):
			ax=fig0.add_axes([0.05+0.2*i,0.05+0.3*j,0.2,0.3])
			ax.tick_params(pad=4)
			ax2=fig2.add_axes([0.05+0.2*i,0.05+0.3*j,0.2,0.3])
			ax2.tick_params(pad=4)
			ax3=fig3.add_axes([0.05+0.2*i,0.05+0.3*j,0.2,0.3])
			ax3.tick_params(pad=4)
			print np.mean(np.clip(N[I[i]*6+J[j]],1,1e12))
			#print np.mean(np.log10(FF0[(I[i]+3)*6+J[j]]*1./FF0[I[i]*6+J[j]]))
			#print np.mean(np.log10(np.clip(N[(I[i]+3)*6+J[j]]*1./Nff[I[i]+3],1e-6,1e12)*1./np.clip(N[I[i]*6+J[j]]\
			#*1./Nff[I[i]],1e-6,1e12)))
			X=np.log10(FF0[(I[i]+2)*6+J[j]]*1./FF0[I[i]*6+J[j]])\
			*1./np.clip(np.log10(np.clip(N[(I[i]+2)*6+J[j]]*1./Nff[I[i]+2],1e-6,1e12)*1./np.clip(N[I[i]*6+J[j]]\
			*1./Nff[I[i]],1e-6,1e12)),1e-3,1e12)
			FF2=np.reshape(FF0[I[i]*6+J[j]].copy(),(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]]))
			FF2[3:,3:]=np.repeat(np.repeat(np.reshape(FF0[I[i]*6+J2[j]],(nB[I[i]*6+J2[j]],nB[I[i]*6+J2[j]]))\
			,2**(J2[j]-J[j]),axis=0),2**(J2[j]-J[j]),axis=1)
			N2=np.reshape(N[I[i]*6+J[j]].copy(),(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]]))
			N2[3:,3:]=np.repeat(np.repeat(np.reshape(N[I[i]*6+J2[j]],(nB[I[i]*6+J2[j]],nB[I[i]*6+J2[j]]))\
			,2**(J2[j]-J[j]),axis=0),2**(J2[j]-J[j]),axis=1)
			X2=np.log10(FF2.flatten()*1./FF0[I[i]*6+J[j]])\
			*1./np.clip(np.log10(np.clip(N2.flatten(),1,1e12)*1./np.clip(N[I[i]*6+J[j]],1,1e12)),1e-6,1e12)
			ax.imshow(np.reshape(np.ma.array(X,mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]])),vmin=-0.5, vmax=1.\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			print np.mean(X),np.mean(X2), np.mean(N[I[i]*6+J[j]]),np.mean(N2), np.mean(FF0[I[i]*6+J[j]]), np.mean(FF2)
			ax.set_xticks(np.array([]))
			ax.set_yticks(np.array([]))
			if i==0:
				ax.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax.set_xlabel(Tlabel[i]+' s')
			ax2.imshow(np.reshape(np.ma.array(X2,mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]])),vmin=-0.2, vmax=1.3\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			ax2.set_xticks(np.array([]))
			ax2.set_yticks(np.array([]))
			if i==0:
				ax2.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax2.set_xlabel(Tlabel[i]+' s')
			ax3.imshow(np.reshape(np.ma.array(X-X2,mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]])),vmin=-1.3, vmax=0.7\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			ax3.set_xticks(np.array([]))
			ax3.set_yticks(np.array([]))
			if i==0:
				ax3.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax3.set_xlabel(Tlabel[i]+' s')
	ax0=fig0.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.5),10**(1.)))
	cbar0.set_label(r'Fano Factor slope (T)',fontsize='x-small')
	ax0=fig2.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.2),10**(1.3)))
	cbar0.set_label(r'Fano Factor slope (area)',fontsize='x-small')
	ax0=fig3.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-1.3),10**(0.7)))
	cbar0.set_label(r'Fano Factor slope diff',fontsize='x-small')
	pylab.figure(1)
	pylab.savefig(FileName+'_ldiffT.png')
	pylab.close()
	pylab.figure(2)
	pylab.savefig(FileName+'_ldiffA.png')
	pylab.close()
	pylab.figure(3)
	pylab.savefig(FileName+'_ldiffT_A.png')
	pylab.close()
	return

def FFxtPlot4(HdfFile, FileName):
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=g['tMax'].value
	nFrames=g['nFrames'].value
	g.close()
	g=h5py.File(FileName+'.hdf','r+')
	Nff=g['Nff'].value
	Sq=[]
	N=[]
	for i in range(14):
		for j in range(6):
			Sq.append(g['Sq%i_%i'%(i,j)].value)
			N.append(g['N%i_%i'%(i,j)].value)
	g.close()
	Res=8
	nBins=64*Res
	Area=4.**np.arange(-1,5)# in number of electrodes (j)
	T=20./Sampling*2.**np.arange(14)
	Tlabel=np.array(['%1.4f'%T[8],'%1.3f'%T[9],'%1.2f'%T[10],'%2.1f'%T[11]])
	Alabel=np.array(['0.25','1','4'])
	NbyA=[]
	FF0=[]
	nB=[]
	TbyA=np.zeros(6*14)
	for i in range(14):
		#nB.append(nBins-3)
		for j in range(6):
			#if j>0:
			nB.append(nBins/2**j-3)
			#NbyA.append((3*np.argsort(np.argsort(N[i*6+j])))/np.size(N[i*6+j]))#classify into 3 groups
			TbyA[i*6+j]=(i-j+5)#colourscale[0:24)#T[i]/A[j]#has discrete values...use log and colourscale?
			FF0.append(np.clip((Sq[i*6+j]-N[i*6+j]**2*1./np.clip(Nff[i],1,1e12))*1./np.clip(N[i*6+j],1,1e12),1e-6,1e12))
			#NbyA.append((3*np.argsort(np.argsort(FF0[i*6+j])))/np.size(FF0[i*6+j]))#classify into 3 groups
	fig0=pylab.figure(1,figsize=(9,6))
	fig2=pylab.figure(2,figsize=(9,6))
	fig3=pylab.figure(3,figsize=(9,6))
	#pylab.subplots_adjust(left=0.12, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3)
	I=np.array([8,11,10,11])
	J=np.array([0,1,2])
	J2=np.array([1,2,3])
	#print N
	#for i in range(30):
	#	print FF0[i].min(), FF0[i].max(), NbyA[i].mean(), Sq[i].mean(), N[i].mean()
	for i in range(2):
		for j in range(3):
			ax=fig0.add_axes([0.05+0.4*i,0.05+0.3*j,0.2,0.3])
			ax.tick_params(pad=4)
			ax2=fig2.add_axes([0.05+0.4*i,0.05+0.3*j,0.2,0.3])
			ax2.tick_params(pad=4)
			ax2a=fig2.add_axes([0.05+0.4*i+0.2,0.05+0.3*j,0.2,0.3])
			ax2a.tick_params(pad=4)
			ax3=fig3.add_axes([0.05+0.4*i,0.05+0.3*j,0.2,0.3])
			ax3.tick_params(pad=4)
			print np.mean(np.clip(N[I[i]*6+J[j]],1,1e12))
			#print np.mean(np.log10(FF0[(I[i]+3)*6+J[j]]*1./FF0[I[i]*6+J[j]]))
			#print np.mean(np.log10(np.clip(N[(I[i]+3)*6+J[j]]*1./Nff[I[i]+3],1e-6,1e12)*1./np.clip(N[I[i]*6+J[j]]\
			#*1./Nff[I[i]],1e-6,1e12)))
			X=np.log10(FF0[(I[i]+2)*6+J[j]]*1./FF0[I[i]*6+J[j]])\
			*1./np.clip(np.log10(np.clip(N[(I[i]+2)*6+J[j]]*1./Nff[I[i]+2],1e-6,1e12)*1./np.clip(N[I[i]*6+J[j]]\
			*1./Nff[I[i]],1e-6,1e12)),1e-3,1e12)
			FF2=np.reshape(FF0[I[i]*6+J[j]].copy(),(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]]))
			FF2[3:,3:]=np.repeat(np.repeat(np.reshape(FF0[I[i]*6+J2[j]],(nB[I[i]*6+J2[j]],nB[I[i]*6+J2[j]]))\
			,2**(J2[j]-J[j]),axis=0),2**(J2[j]-J[j]),axis=1)
			N2=np.reshape(N[I[i]*6+J[j]].copy(),(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]]))
			N2[3:,3:]=np.repeat(np.repeat(np.reshape(N[I[i]*6+J2[j]],(nB[I[i]*6+J2[j]],nB[I[i]*6+J2[j]]))\
			,2**(J2[j]-J[j]),axis=0),2**(J2[j]-J[j]),axis=1)
			X2=np.log10(FF2.flatten()*1./FF0[I[i]*6+J[j]])\
			*1./np.clip(np.log10(np.clip(N2.flatten(),1,1e12)*1./np.clip(N[I[i]*6+J[j]],1,1e12)),1e-6,1e12)
			X2r=np.reshape(X2,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]]))
			X2o=np.zeros((nB[I[i]*6+J[j]]-10,nB[I[i]*6+J[j]]-10))
			X2v=np.cumsum(X2r,axis=0)
			W=np.array([2/16.,-16/20.**2,16/28.**2,-16/36.**2])
			for ii in (4,6,8,10):
				X2h=np.cumsum(X2v[ii:,:]-X2v[:-ii,:],axis=1)
				print X2h.shape
				X2o+=W[(ii-4)/2]*(X2h[:,ii:]-X2h[:,:-ii])\
				[(10-ii)/2:nB[I[i]*6+J[j]]-5-ii/2,(10-ii)/2:nB[I[i]*6+J[j]]-5-ii/2]
			ax.imshow(np.reshape(np.ma.array(X,mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]])),vmin=-0.5, vmax=1.\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			print np.mean(X),np.mean(X2), np.mean(N[I[i]*6+J[j]]),np.mean(N2), np.mean(FF0[I[i]*6+J[j]]), np.mean(FF2)
			ax.set_xticks(np.array([]))
			ax.set_yticks(np.array([]))
			if i==0:
				ax.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax.set_xlabel(Tlabel[i]+' s')
			ax2.imshow(np.reshape(np.ma.array(X2,mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]])),vmin=-0.2, vmax=1.3\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			ax2a.imshow(X2o,vmin=-0.2, vmax=1.3\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			ax2.set_xticks(np.array([]))
			ax2.set_yticks(np.array([]))
			ax2a.set_xticks(np.array([]))
			ax2a.set_yticks(np.array([]))
			if i==0:
				ax2.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax2.set_xlabel(Tlabel[i]+' s')
			ax3.imshow(np.reshape(np.ma.array(X-X2,mask=(N[I[i]*6+J[j]]<10))\
			,(nB[I[i]*6+J[j]],nB[I[i]*6+J[j]])),vmin=-1.3, vmax=0.7\
			,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
			,origin='lower',extent=(0,nB[I[i]*6+J[j]],0,nB[I[i]*6+J[j]]))
			ax3.set_xticks(np.array([]))
			ax3.set_yticks(np.array([]))
			if i==0:
				ax3.set_ylabel(Alabel[j]+' electrodes')
			if j==0:
				ax3.set_xlabel(Tlabel[i]+' s')
	ax0=fig0.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.5),10**(1.)))
	cbar0.set_label(r'Fano Factor slope (T)',fontsize='x-small')
	ax0=fig2.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.2),10**(1.3)))
	cbar0.set_label(r'Fano Factor slope (area)',fontsize='x-small')
	ax0=fig3.add_axes([0.9, 0.1, 0.025, 0.8])
	ax0.tick_params(pad=4)
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-1.3),10**(0.7)))
	cbar0.set_label(r'Fano Factor slope diff',fontsize='x-small')
	pylab.figure(1)
	#pylab.savefig(FileName+'_ldiffT.png')
	pylab.close()
	pylab.figure(2)
	pylab.savefig(FileName+'_avgdiffA.png')
	pylab.close()
	pylab.figure(3)
	#pylab.savefig(FileName+'_ldiffT_A.png')
	pylab.close()
	return

def FFxtPlot5(HdfFile, FileName):
	#f=h5py.File(NoisyChFile,'r+')
	#NoisyAreas=f['NoisyAreas'].value
	#f.close()
	#nBins0=int(np.sqrt(NoisyAreas.shape[0]))
	g=h5py.File(HdfFile,'r+')
	Sampling=g['Sampling'].value
	tMax=int(g['tMax'].value)
	nFrames=g['nFrames'].value
	SingleSpk=np.array(g['RepolarizingSpikes'].value,dtype=bool)
	Loc=(g['Locations'].value)[SingleSpk]
	Amp=g['Amplitudes'].value[SingleSpk]
	Times=np.array((g['Times'].value)[SingleSpk],dtype=int)
	#Loc=Loc[Amp>3.,:]
	#Times=Times[Amp>3.]
	#Times=np.concatenate((Times,np.array([tMax*Sampling+1])))
	I=np.concatenate((np.array([0]),np.cumsum(np.histogram(Times,bins=np.arange(0,int(tMax*Sampling)+20,20))[0])))
	#NNoise=g['CorrelationAnalysis/Noise'].value
	#nBins0=g['CorrelationAnalysis/Parameter/Resolution'].value
	Res=2
	Res0=4
	dt=2
	nRep=4
	nBins=64*Res
	Sm=np.zeros((nBins*Res0,nBins*Res0))
	Ssq=np.zeros((nBins*Res0,nBins*Res0))
	Sm2=np.zeros((nBins*Res0,nBins*Res0))
	Ssq2=np.zeros((nBins*Res0,nBins*Res0))
	for k in range(Res0):
		for kk in range(Res0):
			print kk
			Spikes=np.clip(np.array(Res*(Loc[:,0]-(k+0.5)/Res0+1./2),dtype=int),0,nBins-1)*nBins\
			+np.clip(np.array(Res*(Loc[:,1]-(kk+0.5)/Res0+1./2),dtype=int),0,nBins-1)
			for i in range(nRep):
				S=np.histogram2d(Spikes,Times,bins=(np.arange(nBins**2+1),np.linspace(0+i*tMax/dt/nRep\
				,tMax*Sampling+(nRep-i)*tMax/dt/nRep,tMax/dt+1)))[0]
				S2=np.reshape(np.repeat(np.repeat(np.sum(np.sum(np.reshape(S,(nBins/2,2,nBins/2,2,-1))\
				,axis=3),axis=1),2,axis=0),2,axis=1),(nBins**2,-1))#misaligned!!!(shift by Res0/2)
				Sm[k::Res0,kk::Res0]+=np.reshape(np.sum(S,axis=1),(nBins,nBins))
				Ssq[k::Res0,kk::Res0]+=np.reshape(np.sum(S**2,axis=1),(nBins,nBins))
				Sm2[k::Res0,kk::Res0]+=np.reshape(np.sum(S2,axis=1),(nBins,nBins))
				Ssq2[k::Res0,kk::Res0]+=np.reshape(np.sum(S2**2,axis=1),(nBins,nBins))
				#print S.mean(),S2.mean(), Sm.mean(), Sm2.mean(), Ssq.mean(),Ssq2.mean()
	FF0=np.clip((Ssq-Sm**2*1./(nRep*(tMax/dt)))*1./np.clip(Sm,1,1e12),1e-3,1e12)#1
	FF2=np.clip((Ssq2-Sm2**2*1./(nRep*(tMax/dt)))*1./np.clip(Sm2,1,1e12),1e-3,1e12)#2
	X2=np.log10((FF2)*1./(FF0))\
	*1./np.clip(np.log10(np.clip(Sm2,2,1e12)\
	*1./np.clip(Sm,1,1e12)),1e-3,1e12)#3, 2*Res0 shorter
	print X2.mean(), X2.var()
	fig0=pylab.figure(1,figsize=(7,4))
	ax3=fig0.add_axes([0.0,0.1,0.2,0.35])
	ax0=fig0.add_axes([0.0,0.1+0.45,0.2,0.35])
	ax4=fig0.add_axes([0.0+0.33,0.1,0.2,0.35])
	ax1=fig0.add_axes([0.0+0.33,0.1+0.45,0.2,0.35])
	ax5=fig0.add_axes([0.0+0.66,0.1,0.2,0.35])
	ax2=fig0.add_axes([0.0+0.66,0.1+0.45,0.2,0.35])
	ax3a=fig0.add_axes([0.2,0.1,0.025,0.35])
	ax0a=fig0.add_axes([0.2,0.1+0.45,0.025,0.35])
	ax4a=fig0.add_axes([0.3+0.23,0.1,0.025,0.35])
	ax1a=fig0.add_axes([0.3+0.23,0.1+0.45,0.025,0.35])
	ax5a=fig0.add_axes([0.3+0.56,0.1,0.025,0.35])
	ax2a=fig0.add_axes([0.3+0.56,0.1+0.45,0.025,0.35])
	ax0.tick_params(pad=4)
	ax1.tick_params(pad=4)
	ax2.tick_params(pad=4)
	ax3.tick_params(pad=4)
	ax4.tick_params(pad=4)
	ax5.tick_params(pad=4)
	ax0.imshow(np.log10(FF0),vmin=-0.5, vmax=3.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
	,origin='lower',extent=(0,64,0,64))
	ax0.set_xticks(np.array([0,64]))
	ax0.set_yticks(np.array([0,64]))
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax0a,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.5),10**(3.)))
	cbar0.set_label(r'Fano Factor (0.25el)',fontsize='x-small')
	ax1.imshow(np.log10(FF2),vmin=-0.5, vmax=3.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
	,origin='lower',extent=(0,64,0,64))
	ax1.set_xticks(np.array([0,64]))
	ax1.set_yticks(np.array([0,64]))
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax1a,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.5),10**(3.)))
	cbar0.set_label(r'Fano Factor (1el)',fontsize='x-small')
	ax2.imshow(X2,vmin=-0.5, vmax=2.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
	,origin='lower',extent=(0,64,0,64))
	ax2.set_xticks(np.array([0,64]))
	ax2.set_yticks(np.array([0,64]))
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax2a,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.5),10**(2.)))
	cbar0.set_label(r'Fano Factor ratio (area)',fontsize='x-small')
	X2v=np.cumsum(np.concatenate((np.zeros(nBins*Res0)[None,:],X2),axis=0),axis=0)
	MaxXo=np.zeros((nBins*Res0-1,nBins*Res0-1))
	MinXo=np.zeros((nBins*Res0-1,nBins*Res0-1))
	for ii in (4,6):
		X2o=np.zeros((nBins*Res0-ii-1,nBins*Res0-ii-1))
		X2h=np.cumsum(np.concatenate((np.zeros(nBins*Res0-ii-1)[:,None],X2v[1+ii:-1,:]-X2v[1:-ii-1,:]),axis=1),axis=1)
		X2h2=np.cumsum(np.concatenate((np.zeros(nBins*Res0-ii-1)[:,None],X2v[ii+2:,:]-X2v[:-ii-2,:]),axis=1),axis=1)
		print X2h.shape
		X2o+=1./ii**2+1./(ii+2.)**2*(X2h[:,1+ii:-1]-X2h[:,1:-ii-1])
		X2o+=1./(ii+2.)**2*(X2h2[:,ii+2:]-X2h2[:,:-ii-2])
		MaxXo[ii/2:-ii/2,ii/2:-ii/2]=0.5*np.abs(MaxXo[ii/2:-ii/2,ii/2:-ii/2]-X2o)\
		+(MaxXo[ii/2:-ii/2,ii/2:-ii/2]+X2o)/2.
		MinXo[ii/2:-ii/2,ii/2:-ii/2]=-0.5*np.abs(MinXo[ii/2:-ii/2,ii/2:-ii/2]-X2o)\
		+(MinXo[ii/2:-ii/2,ii/2:-ii/2]+X2o)/2.
	ax3.imshow(MaxXo,vmin=-0., vmax=2.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
	,origin='lower',extent=(0,64.-(1./Res0)/Res,0,64.-(1./Res0)/Res))
	ax3.set_xticks(np.array([0,63]))
	ax3.set_yticks(np.array([0,63]))
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax3a,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.),10**(2.)))
	cbar0.set_label(r'Fano Factor local max.',fontsize='x-small')
	ax4.imshow(MinXo,vmin=-0.5, vmax=0.5\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
	,origin='lower',extent=(0,64.-(1./Res0)/Res,0,64.-(1./Res0)/Res))
	ax4.set_xticks(np.array([0,63]))
	ax4.set_yticks(np.array([0,63]))
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax4a,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.5),10**(0.5)))
	cbar0.set_label(r'Fano Factor local min.',fontsize='x-small')
	ax5.imshow(MaxXo-MinXo,vmin=-0., vmax=2.\
	,cmap=pylab.cm.spectral, aspect='equal',interpolation='none'\
	,origin='lower',extent=(0,64.-(1./Res0)/Res,0,64.-(1./Res0)/Res))
	ax5.set_xticks(np.array([0,63]))
	ax5.set_yticks(np.array([0,63]))
	cbar0 =mplt.colorbar.ColorbarBase(ax=ax5a,cmap=pylab.cm.spectral, norm=mplt.colors.LogNorm(10**(-0.),10**(2.)))
	cbar0.set_label(r'Fano Factor local variability',fontsize='x-small')
	pylab.figure(1)
	pylab.savefig(FileName+'_locMax_2s.png')
	pylab.close()
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
	fig2=pylab.figure(2,figsize=(7.,3.5))
	fig3=pylab.figure(3,figsize=(7.,3.2))
	#pylab.subplots_adjust(left=0.12, bottom=0.1, right=0.95, top=0.9, wspace=0.2, hspace=0.3)
	ax0=fig2.add_axes([0.05, 0.1, 0.38, 0.76])
	ax3=fig3.add_axes([0.05, 0.08, 0.35, 0.78])
	#ax2=fig3.add_axes([0.63, 0.54, 0.25, 0.38])
	ax4=fig3.add_axes([0.45, 0.08, 0.35, 0.78])
	ax1=fig2.add_axes([0.48, 0.1, 0.38, 0.76])
	ax5=fig3.add_axes([0.81, 0.08, 0.18, 0.45])
	ax6=fig2.add_axes([0.87, 0.1, 0.02, 0.76])
	#ax7=fig3.add_axes([0.92, 0.05, 0.015, 0.38])
	ax0.tick_params(pad=8)
	ax1.tick_params(pad=8)
	#ax2 .tick_params(pad=8)
	ax3.tick_params(pad=8)
	ax4.tick_params(pad=8)
	#ax5.tick_params(pad=8)
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

def TempBiasPlot(HdfFile, FileName):
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
	TempBias=g['TempBias'].value#should be 0.5 on average
	print np.min(TempBias), np.max(TempBias), np.mean(TempBias)
	nBins=Res*64
	#convert to colour
	a=np.array([[1,0,0],[0.5,0.5,0],[0,1,0],[0,0.5,0.5],[0,0,1]])
	cDloc1=np.dot(Dloc1,a)
	cDloc1=cDloc1*1./np.clip(np.sum(cDloc1,axis=2),1e-6,1e12)[:,:,None]
	#need to compute something related to significance
	cDavg=np.clip(np.dot((Dloc1-Davg)*np.sqrt(Nloc1)[:,:,None]*Ns/3.,a)+0.5,0,1)
	#modify such that bias reflects brightness
	cDavg*=(1.-2.*np.abs(TempBias-0.5))[:,:,None]
	cDavg+=(np.abs(TempBias-0.5))[:,:,None]
	cDavg+=(TempBias-0.5)[:,:,None]
	b=np.array([[1,0,0],[0.5,0.5,0],[0,1,0],[0,0.5,0.5],[0,0,1],[0.5,0,0.5]])
	cAloc1=np.dot(Aloc1,b)
	cAloc1=cAloc1*1./np.clip(np.sum(cAloc1,axis=2),1e-6,1e12)[:,:,None]
	cAavg=np.clip(np.dot((Aloc1-Aavg)*np.sqrt(Nloc1)[:,:,None]*Ns/3.,b)+0.5,0,1)
	#modify such that bias reflects brightness
	cAavg*=(1.-2.*np.abs(TempBias-0.5))[:,:,None]
	cAavg+=(np.abs(TempBias-0.5))[:,:,None]
	cAavg+=(TempBias-0.5)[:,:,None]
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
	ax0.tick_params(pad=8)
	ax1.tick_params(pad=8)
	#ax2 .tick_params(pad=8)
	ax3.tick_params(pad=8)
	ax4.tick_params(pad=8)
	#ax5.tick_params(pad=8)
	fig2.text(0.01,0.91,'A',fontsize='xx-large')
	fig2.text(0.41,0.91,'B',fontsize='xx-large')
	fig3.text(0.01,0.91,'A',fontsize='xx-large')
	fig3.text(0.41,0.91,'B',fontsize='xx-large')
	#fig3.text(0.31,0.45,'E',fontsize='xx-large')
	#fig3.text(0.61,0.45,'F',fontsize='xx-large')
	ax0.imshow(TempBias*4.-1.5,vmin=0., vmax=1., cmap=pylab.cm.spectral ,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,nBins-Ns+1,0,nBins-Ns+1))
	ax0.set_xlim(-2,nBins-2)
	ax0.set_ylim(-2,nBins-2)
	ax0.set_xticks(np.array([0,nBins/2,nBins])-2)
	ax0.set_xticklabels(np.array([0,32,64]))
	ax0.set_yticks(np.array([0,nBins/2,nBins])-2)
	ax0.set_yticklabels(np.array([0,32,64]))
	ax3.imshow(TempBias,vmin=0., vmax=1., cmap=pylab.cm.spectral ,aspect='equal'\
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


def FRBiasPlot(HdfFile, FileName):
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
	TempBias=g['FRBias'].value#should be 0.5 on average
	print np.min(TempBias), np.max(TempBias), np.mean(TempBias)
	nBins=Res*64
	#convert to colour
	a=np.array([[1,0,0],[0.5,0.5,0],[0,1,0],[0,0.5,0.5],[0,0,1]])
	cDloc1=np.dot(Dloc1,a)
	cDloc1=cDloc1*1./np.clip(np.sum(cDloc1,axis=2),1e-6,1e12)[:,:,None]
	#need to compute something related to significance
	cDavg=np.clip(np.dot((Dloc1-Davg)*np.sqrt(Nloc1)[:,:,None]*Ns/3.,a)+0.5,0,1)
	#modify such that bias reflects brightness
	cDavg*=(1.-2.*np.abs(TempBias-0.5))[:,:,None]
	cDavg+=(np.abs(TempBias-0.5))[:,:,None]
	cDavg+=(TempBias-0.5)[:,:,None]
	b=np.array([[1,0,0],[0.5,0.5,0],[0,1,0],[0,0.5,0.5],[0,0,1],[0.5,0,0.5]])
	cAloc1=np.dot(Aloc1,b)
	cAloc1=cAloc1*1./np.clip(np.sum(cAloc1,axis=2),1e-6,1e12)[:,:,None]
	cAavg=np.clip(np.dot((Aloc1-Aavg)*np.sqrt(Nloc1)[:,:,None]*Ns/3.,b)+0.5,0,1)
	#modify such that bias reflects brightness
	cAavg*=(1.-2.*np.abs(TempBias-0.5))[:,:,None]
	cAavg+=(np.abs(TempBias-0.5))[:,:,None]
	cAavg+=(TempBias-0.5)[:,:,None]
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
	ax0.tick_params(pad=8)
	ax1.tick_params(pad=8)
	#ax2 .tick_params(pad=8)
	ax3.tick_params(pad=8)
	ax4.tick_params(pad=8)
	#ax5.tick_params(pad=8)
	fig2.text(0.01,0.91,'A',fontsize='xx-large')
	fig2.text(0.41,0.91,'B',fontsize='xx-large')
	fig3.text(0.01,0.91,'A',fontsize='xx-large')
	fig3.text(0.41,0.91,'B',fontsize='xx-large')
	#fig3.text(0.31,0.45,'E',fontsize='xx-large')
	#fig3.text(0.61,0.45,'F',fontsize='xx-large')
	ax0.imshow(TempBias*4.-1.5,vmin=0., vmax=1., cmap=pylab.cm.spectral ,aspect='equal'\
	,interpolation='none',origin='lower',extent=(0,nBins-Ns+1,0,nBins-Ns+1))
	ax0.set_xlim(-2,nBins-2)
	ax0.set_ylim(-2,nBins-2)
	ax0.set_xticks(np.array([0,nBins/2,nBins])-2)
	ax0.set_xticklabels(np.array([0,32,64]))
	ax0.set_yticks(np.array([0,nBins/2,nBins])-2)
	ax0.set_yticklabels(np.array([0,32,64]))
	ax3.imshow(TempBias,vmin=0., vmax=1., cmap=pylab.cm.spectral ,aspect='equal'\
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

def TransitionMatricesPlot(HdfFile1, HdfFile2, FileName\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(100.*np.exp(-cQ1)))*100./len(cQ1)
	A2=np.argsort(np.argsort(100.*np.exp(-cQ2)))*100./len(cQ2)
	L1=len(A1)
	L2=len(A2)
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	#H1r=np.reshape(H1,(-1,100))
	#H2r=np.reshape(H2,(-1,100))
	#histogram of spike count
	HN1=np.histogram2d(Loc1[:,0],Loc1[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN2=np.histogram2d(Loc2[:,0],Loc2[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN1=HN1[1:,1:]+HN1[:-1,1:]+HN1[:-1,:-1]+HN1[1:,:-1]
	HN2=HN2[1:,1:]+HN2[:-1,1:]+HN2[:-1,:-1]+HN2[1:,:-1]
	HN1=HN1.flatten()
	HN2=HN2.flatten()
	#threshold spike count
	#NInd=(HN1+HN2)>=100
	#NInd=np.argsort(HN1+HN2)[-500:]
	#Hsum=np.clip(np.sum((H2r-H1r)*np.arange(100)[None,:],axis=1)*1./np.clip((HN2+HN1),1,1e12),-100,100)
	##Hloc=np.histogram(Hsum,bins=np.arange(-100,101))[0]
	#Histogram of average and single location activity change
	#print Hsum.shape
	#Q=np.histogram2d((Hsum[:,None]*np.ones(100)[None,:]).flatten()\
	#,(np.ones(Hsum.shape[0])[:,None]*np.arange(100)[None,:]).flatten()\
	#,bins=(np.arange(-100,101),np.arange(101)),weights=(H2r-H1r).flatten())[0]
	#QN=np.histogram(Hsum,bins=(np.arange(-100,101)),weights=np.sqrt(np.sum(H2r+H1r,axis=1)))[0]#
	#need to figure out where the largest difference is (in Ex.) take 10% window (sliding, periodic boundary)
	#use absolute difference (^1/3?)
	#cannot plot that...clustering? -- distances? - correlations... or firing rate cutoff take highest 500?
	#S=((H2r-H1r)*1./np.clip((np.sqrt(np.sum(H2r+H1r,axis=1))),1,1e12)[:,None])[NInd,:][np.argsort(Hsum[NInd]),:]
	
	
	Rgbmap=np.zeros((100,3))
	Rgbmap[17:83,1]=np.bartlett(66)
	Rgbmap[:49,0]=Rgbmap[34:83,1]
	Rgbmap[-17:,0]=Rgbmap[17:34,1]
	Rgbmap[-49:,2]=Rgbmap[17:66,1]
	Rgbmap[:17,2]=Rgbmap[66:83,1]
	#find difference maxima
	D=np.cumsum(np.concatenate(((H2-H1)[:,:,-6:],(H2-H1),(H2-H1)[:,:,:5]),axis=2),axis=2)
	A=np.cumsum(np.concatenate(((H2+H1)[:,:,-6:],(H2+H1),(H2+H1)[:,:,:5]),axis=2),axis=2)
	#print D.shape, H1.shape, H2.shape
	dD=D[:,:,11:]-D[:,:,:-11]
	dA=np.clip(np.sqrt(A[:,:,11:]-A[:,:,:-11]),1,1e12)
	a=np.argmax(np.abs(dD/dA),axis=2).flatten()
	#b=np.zeros(a.shape)
	b=np.clip(((dD/dA).flatten()[a+np.arange(dD.shape[0]*dD.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	#convert to color
	X=np.reshape(Rgbmap[a,:],((yMax-yMin)*Res,-1,3))
	Y=(X-0.5)*np.reshape(2*b,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	#(H2-H1)/(H2+H1)
	
	W=np.clip(np.dot(np.reshape((dD/dA),(-1,100))-0.5,Rgbmap-0.5)/150.+0.5,0,1)
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (2./12., 1., 1.), (6./12., 0.0, 0.0),\
	(10./12., 0.0, 0.0), (1.0, 0.75, 0.75)),\
	'green': ((0.0,0.,0.),(2./12.,0.,0.),(6./12.,1.,1.),(10./12.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(2./12.,0.,0.),(6./12.,0.,0.),(10./12.,1.,1.),(1.0, 0.25, 0.25))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (2./12., 0.0, 0.0), (6./12., 1., 1.),\
	(10./12., 1., 1.), (1.0, 0.25, 0.25)),\
	'green': ((0.0,1.,1.),(1./12.,1.,1.),(6./12.,0.0,0.0),(10./12.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(2./12.,1.,1.),(6./12.,1.,1.),(10./12.,0.0,0.0),(1.0, 0.75, 0.75))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	
	
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (decreasing)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	#cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (increasing)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	#ax.imshow((Q*1./(QN+(QN==0))[:,None]).transpose()\
	#ax.imshow(S.transpose()\
	#,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	#,interpolation='none',origin='lower',extent=(-100,100,0,100))
	#ax.set_ylim(0,100)
	#ax.set_xlim(-100,100)
	#ax.set_xticks(np.array([-100,0,100]))
	#ax.set_yticks(np.array([0,50,100]))
	ax.imshow(np.reshape(W,((yMax-yMin)*Res,-1,3)), aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(Y\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	#pylab.figure()
	#pylab.imshow(np.reshape(Hsum,((yMax-yMin)*Res,-1))\
	#,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	#pylab.colorbar()
	#pylab.show()
	return

def TransitionMatricesPlotDev(HdfFile1, HdfFile2, FileName, FileName2\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(100.*np.exp(-cQ1)))*100./len(cQ1)
	A2=np.argsort(np.argsort(100.*np.exp(-cQ2)))*100./len(cQ2)
	L1=len(A1)
	L2=len(A2)
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	#H1r=np.reshape(H1,(-1,100))
	#H2r=np.reshape(H2,(-1,100))
	#histogram of spike count
	HN1=np.histogram2d(Loc1[:,0],Loc1[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN2=np.histogram2d(Loc2[:,0],Loc2[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN1=HN1[1:,1:]+HN1[:-1,1:]+HN1[:-1,:-1]+HN1[1:,:-1]
	HN2=HN2[1:,1:]+HN2[:-1,1:]+HN2[:-1,:-1]+HN2[1:,:-1]
	HN1=HN1.flatten()
	HN2=HN2.flatten()
	#threshold spike count
	#NInd=(HN1+HN2)>=100
	#NInd=np.argsort(HN1+HN2)[-500:]
	#Hsum=np.clip(np.sum((H2r-H1r)*np.arange(100)[None,:],axis=1)*1./np.clip((HN2+HN1),1,1e12),-100,100)
	##Hloc=np.histogram(Hsum,bins=np.arange(-100,101))[0]
	#Histogram of average and single location activity change
	#print Hsum.shape
	#Q=np.histogram2d((Hsum[:,None]*np.ones(100)[None,:]).flatten()\
	#,(np.ones(Hsum.shape[0])[:,None]*np.arange(100)[None,:]).flatten()\
	#,bins=(np.arange(-100,101),np.arange(101)),weights=(H2r-H1r).flatten())[0]
	#QN=np.histogram(Hsum,bins=(np.arange(-100,101)),weights=np.sqrt(np.sum(H2r+H1r,axis=1)))[0]#
	#need to figure out where the largest difference is (in Ex.) take 10% window (sliding, periodic boundary)
	#use absolute difference (^1/3?)
	#cannot plot that...clustering? -- distances? - correlations... or firing rate cutoff take highest 500?
	#S=((H2r-H1r)*1./np.clip((np.sqrt(np.sum(H2r+H1r,axis=1))),1,1e12)[:,None])[NInd,:][np.argsort(Hsum[NInd]),:]
	
	
	Rgbmap=np.zeros((100,3))
	Rgbmap[10:90,1]=np.bartlett(80)
	Rgbmap[10:50,0]=np.bartlett(80)[40:]
	Rgbmap[:10,0]=np.bartlett(80)[30:40]
	Rgbmap[-10:,0]=np.bartlett(80)[:10]
	Rgbmap[50:-10,2]=np.bartlett(80)[:40]
	Rgbmap[-10:,2]=np.bartlett(80)[40:50]
	Rgbmap[:10,2]=np.bartlett(80)[70:]
	#find difference maxima
	D=np.cumsum(np.concatenate(((H2-H1)[:,:,:6][:,:,::-1],(H2-H1),(H2-H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#need to subtract a uniform distribution here...
	D1=np.cumsum(np.concatenate(((H1-1.5*np.mean(H1,axis=2)[:,:,None])[:,:,:6][:,:,::-1]\
	,(H1-1.5*np.mean(H1,axis=2)[:,:,None]),(H1-1.5*np.mean(H1,axis=2)[:,:,None])[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	A=np.cumsum(np.concatenate(((H2+H1)[:,:,:6][:,:,::-1],(H2+H1),(H2+H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#print D.shape, H1.shape, H2.shape
	dD=D[:,:,11:]-D[:,:,:-11]
	dD1=D1[:,:,11:]-D1[:,:,:-11]
	dA=np.clip(np.sqrt(A[:,:,11:]-A[:,:,:-11]),1,1e12)
	a=np.argmax(np.abs(dD/dA),axis=2).flatten()
	a1=np.argmax(np.abs(dD1),axis=2).flatten()
	#a1=np.argsort(np.argsort(np.abs(dD1),axis=2),axis=2).flatten()\
	#[a+np.arange(dD.shape[0]*dD.shape[1],dtype=int)*100]
	#b=np.zeros(a.shape)
	b=np.clip(((dD/dA).flatten()[a+np.arange(dD.shape[0]*dD.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	b1=np.clip(((np.sign(dD1)*np.sqrt(np.abs(dD1))).flatten()[a1+np.arange(dD1.shape[0]*dD1.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	#convert to color
	X=np.reshape(Rgbmap[a,:],((yMax-yMin)*Res,-1,3))
	X1=np.reshape(Rgbmap[a1,:],((yMax-yMin)*Res,-1,3))
	Y=(X-0.5)*np.reshape(2*b,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	Y1=(X1-0.5)*np.reshape(2*b1,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	#(H2-H1)/(H2+H1)
	
	#W=np.clip(np.dot(np.reshape((dD/dA),(-1,100))-0.5,Rgbmap-0.5)/150.+0.5,0,1)
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (1./10., 1., 1.), (5./10., 0.0, 0.0),\
	(9./10., 0.0, 0.0), (1.0, 0.25, 0.25)),\
	'green': ((0.0,0.,0.),(1./10.,0.,0.),(5./10.,1.,1.),(9./10.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(1./10.,0.,0.),(5./10.,0.,0.),(9./10.,1.,1.),(1.0, 0.75, 0.75))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (1./10., 0.0, 0.0), (5./10., 1., 1.),\
	(9./10., 1., 1.), (1.0, 0.75, 0.75)),\
	'green': ((0.0,1.,1.),(1./10.,1.,1.),(5./10.,0.0,0.0),(9./10.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(1./10.,1.,1.),(5./10.,1.,1.),(9./10.,0.0,0.0),(1.0, 0.25, 0.25))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	
	
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (decreasing)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	#cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (increasing)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	#ax.imshow((Q*1./(QN+(QN==0))[:,None]).transpose()\
	#ax.imshow(S.transpose()\
	#,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	#,interpolation='none',origin='lower',extent=(-100,100,0,100))
	#ax.set_ylim(0,100)
	#ax.set_xlim(-100,100)
	#ax.set_xticks(np.array([-100,0,100]))
	#ax.set_yticks(np.array([0,50,100]))np.reshape(W,((yMax-yMin)*Res,-1,3))
	ax.imshow(Y1, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(Y\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	c=np.histogram2d(np.sign(b)*50+a,np.sign(b1)*50+a1,bins=(np.arange(-50,150),np.arange(-50,150)))[0]
	fig2=pylab.figure(1,figsize=(8,8))
	fig2.text(0.05,0.5,'typical excitability rank',rotation='vertical',ha='center',va='center')#,fontsize='xx-large'
	fig2.text(0.5,0.05,'typical excitability change',ha='center',va='center')#,fontsize='xx-large'
	fig2.text(0.97,0.3,'neg.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	fig2.text(0.97,0.7,'pos.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	fig2.text(0.3,0.97,'neg. change',ha='center',va='center',fontsize='small')
	fig2.text(0.7,0.97,'pos. change',ha='center',va='center',fontsize='small')
	ax4=fig2.add_axes([0.1, 0.1, 0.8, 0.8],axisbg=bgColor)
	ax5=fig2.add_axes([0.1, 0.9, 0.4, 0.03])
	ax6=fig2.add_axes([0.5, 0.9, 0.4, 0.03])
	ax7=fig2.add_axes([0.9, 0.1, 0.03, 0.4])
	ax8=fig2.add_axes([0.9, 0.5, 0.03, 0.4])
	ax4.set_xticks(np.array([]))
	ax4.set_yticks(np.array([]))
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	ax6.set_xticks(np.array([]))
	ax6.set_yticks(np.array([]))
	ax7.set_xticks(np.array([]))
	ax7.set_yticks(np.array([]))
	ax8.set_xticks(np.array([]))
	ax8.set_yticks(np.array([]))
	cbar5 =mplt.colorbar.ColorbarBase(ax=ax5,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar6 =mplt.colorbar.ColorbarBase(ax=ax6,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar7 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar8 =mplt.colorbar.ColorbarBase(ax=ax8,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar5.set_ticks(np.array([]))
	cbar6.set_ticks(np.array([]))
	cbar7.set_ticks(np.array([]))
	cbar8.set_ticks(np.array([]))
	ax4.imshow(np.log10(np.clip(c.transpose(),1,1e3)), aspect='equal', vmin=0., vmax=3., cmap=SpikesCmap\
	, interpolation='none', origin='lower', extent=(0,200,0,200))
	ax4.set_ylim(0,200)
	ax4.set_xlim(0,200)
	#ax4.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	#ax4.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	ax4.plot(np.array([100,100]),np.array([0,200]),'k-')
	ax4.plot(np.array([0,200]),np.array([100,100]),'k-')
	pylab.savefig(FileName2)
	pylab.close()
	#pylab.figure()
	#pylab.imshow(np.reshape(Hsum,((yMax-yMin)*Res,-1))\
	#,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	#pylab.colorbar()
	#pylab.show()
	return

def TransitionMatricesOpt(HdfFile1, HdfFile2, FileName, FileName2\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(100.*np.exp(-cQ1)))*100./len(cQ1)
	A2=np.argsort(np.argsort(100.*np.exp(-cQ2)))*100./len(cQ2)
	L1=len(A1)
	L2=len(A2)
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	D=np.cumsum(np.concatenate(((H2-H1)[:,:,:4][:,:,::-1],(H2-H1),(H2-H1)[:,:,-3:][:,:,::-1]),axis=2),axis=2)
	#D2=np.cumsum(np.concatenate(((H2)[:,:,:6][:,:,::-1],(H2),(H2)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	dD=D[:,:,7:]-D[:,:,:-7]
	N0=dD.shape[0]
	N1=dD.shape[1]
	N2=dD.shape[2]
	DN=np.zeros((N0,N1,N2,4))
	DN[:,:,:,0]=dD.copy()
	for i in range(2,5):
		X=np.cumsum(np.concatenate((np.zeros((1,dD.shape[1],dD.shape[2])),dD),axis=0),axis=0)
		Y=np.cumsum(np.concatenate((np.zeros((dD.shape[0]-i+1,1,dD.shape[2])),X[i:,:,:]-X[:-i,:,:]),axis=1),axis=1)
		DN[(i-1)/2:N0-i/2,(i-1)/2:N1-i/2,:,i-1]=(Y[:,i:,:]-Y[:,:-i,:])*1./i
	RSum=np.zeros((100,4,3,2))
	Rdist=np.zeros((N0,N1,100,4,3))
	#increases in #spikes
	for i in range(4):
		Rdist[:,:,:,i,0]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0]=np.min(np.concatenate(\
					(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:,i][:,:,:,None]),axis=3),axis=3)
		Rdist[:,:,:,i,1]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1]=np.min(np.concatenate(\
					(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:-2,i][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,1:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:,i][:,:,:,None]),axis=3),axis=3)
		Rdist[:,:,:,i,2]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,2:N2-2,i,2]=np.min(np.concatenate(\
					(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,2:N2-2,i,2][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:-4,i][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,3:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:-2,i][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,1:-3,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,4:,i][:,:,:,None]),axis=3),axis=3)
	#need to make that positive
	Rdist=np.clip(Rdist,0,1e12)
	for j in range(1,3):
			Rdist[:,:,:,0,j]=np.min(Rdist[:,:,:,0,:j+1],axis=3)
	for i in range(1,4):
		Rdist[:,:,:,i,0]=np.min(Rdist[:,:,:,i-1:i+1,0],axis=3)
		for j in range(1,3):
			Rdist[:,:,:,i,j]=np.min(np.min(Rdist[:,:,:,:i+1,:j+1],axis=4),axis=3)
	for i in range(4):
		for j in range(3):
			RSum[:,i,j,0]=np.sum(np.sum(Rdist[:,:,:,i,j],axis=1),axis=0)
	RdistN=np.zeros((N0,N1,100,4,3))
	#decreases in #spikes
	for i in range(4):
		RdistN[:,:,:,i,0]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0]=np.max(np.concatenate(\
					(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:,i][:,:,:,None]),axis=3),axis=3)
		RdistN[:,:,:,i,1]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1]=np.max(np.concatenate(\
					(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:-2,i][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,1:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:,i][:,:,:,None]),axis=3),axis=3)
		RdistN[:,:,:,i,2]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,2:N2-2,i,2]=np.max(np.concatenate(\
					(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,2:N2-2,i,2][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:-4,i][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,3:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:-2,i][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,1:-3,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,4:,i][:,:,:,None]),axis=3),axis=3)
	#need to make that negative
	RdistN=np.clip(RdistN,-1e12,0)
	for j in range(1,3):
			RdistN[:,:,:,0,j]=np.max(RdistN[:,:,:,0,:j+1],axis=3)
	for i in range(1,4):
		RdistN[:,:,:,i,0]=np.max(RdistN[:,:,:,i-1:i+1,0],axis=3)
		for j in range(1,3):
			RdistN[:,:,:,i,j]=np.max(np.max(RdistN[:,:,:,:i+1,:j+1],axis=4),axis=3)
	for i in range(4):
		for j in range(3):
			RSum[:,i,j,1]=np.sum(np.sum(RdistN[:,:,:,i,j],axis=1),axis=0)
	Rgbmap=np.zeros((100,3))
	Rgbmap[10:90,1]=np.bartlett(80)
	Rgbmap[10:50,0]=np.bartlett(80)[40:]
	Rgbmap[:10,0]=np.bartlett(80)[30:40]
	Rgbmap[-10:,0]=np.bartlett(80)[:10]
	Rgbmap[50:-10,2]=np.bartlett(80)[:40]
	Rgbmap[-10:,2]=np.bartlett(80)[40:50]
	Rgbmap[:10,2]=np.bartlett(80)[70:]
	#normalization???-- that one seems reasonable...
	A=np.cumsum(np.concatenate(((H2+H1)[:,:,:6][:,:,::-1],(H2+H1),(H2+H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#print D.shape, H1.shape, H2.shape
	dD=Rdist[:,:,:,3,2]
	dD1=RdistN[:,:,:,3,2]
	dA=np.clip(np.sqrt(A[:,:,11:]-A[:,:,:-11]),1,1e12)
	W=np.clip(np.dot(np.reshape((dD/dA),(-1,100))-0.5,Rgbmap-0.5)/50.+0.5,0,1)
	W1=np.clip(np.dot(np.reshape((dD1/dA),(-1,100))-0.5,Rgbmap-0.5)/50.+0.5,0,1)
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (1./10., 1., 1.), (5./10., 0.0, 0.0),\
	(9./10., 0.0, 0.0), (1.0, 0.25, 0.25)),\
	'green': ((0.0,0.,0.),(1./10.,0.,0.),(5./10.,1.,1.),(9./10.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(1./10.,0.,0.),(5./10.,0.,0.),(9./10.,1.,1.),(1.0, 0.75, 0.75))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (1./10., 0.0, 0.0), (5./10., 1., 1.),\
	(9./10., 1., 1.), (1.0, 0.75, 0.75)),\
	'green': ((0.0,1.,1.),(1./10.,1.,1.),(5./10.,0.0,0.0),(9./10.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(1./10.,1.,1.),(5./10.,1.,1.),(9./10.,0.0,0.0),(1.0, 0.25, 0.25))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (decreasing)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (increasing)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	ax.imshow(np.reshape(W,((yMax-yMin)*Res,-1,3)), aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(np.reshape(W1,((yMax-yMin)*Res,-1,3))\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	pylab.figure()
	Clr=np.array(['k','b','g','r'])
	Clr2=np.array(['g','c','y','m'])
	L=np.array(['-','--','-.',':'])
	for i in range(4):
		for j in range(3):
			pylab.plot(np.arange(100),RSum[:,i,j,0],Clr[i]+L[j])
			pylab.plot(np.arange(100),-RSum[:,i,j,1],Clr2[i]+L[j])
	pylab.savefig(FileName2)
	pylab.close()
	return

def TransitionMatricesOpt2(HdfFile1, HdfFile2, FileName, FileName2\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(50.*np.exp(-cQ1)))*50./len(cQ1)
	A2=np.argsort(np.argsort(50.*np.exp(-cQ2)))*50./len(cQ2)
	L1=len(A1)
	L2=len(A2)
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(51)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(51)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	D=np.cumsum(np.concatenate(((H2-H1)[:,:,:3][:,:,::-1],(H2-H1),(H2-H1)[:,:,-2:][:,:,::-1]),axis=2),axis=2)
	dD=D[:,:,5:]-D[:,:,:-5]
	N0=dD.shape[0]
	N1=dD.shape[1]
	N2=dD.shape[2]
	DN=np.zeros((N0,N1,N2,3))
	DN[:,:,:,0]=dD.copy()
	for i in range(2,4):
		X=np.cumsum(np.concatenate((np.zeros((1,dD.shape[1],dD.shape[2])),dD),axis=0),axis=0)
		Y=np.cumsum(np.concatenate((np.zeros((dD.shape[0]-i+1,1,dD.shape[2])),X[i:,:,:]-X[:-i,:,:]),axis=1),axis=1)
		DN[(i-1)/2:N0-i/2,(i-1)/2:N1-i/2,:,i-1]=(Y[:,i:,:]-Y[:,:-i,:])*1./i
	RSum=np.zeros((50,3,2,2))
	Rdist=np.zeros((N0,N1,50,3,2))
	#increases in #spikes
	for i in range(3):
		Rdist[:,:,:,i,0]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0]=np.min(np.concatenate(\
					(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:,i][:,:,:,None]),axis=3),axis=3)
		Rdist[:,:,:,i,1]=DN[:,:,:,i]
		for k in range(i+1):
			for j in range(i+1):
				Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1]=np.min(np.concatenate(\
				(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2,i][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:,i][:,:,:,None]),axis=3),axis=3)
	#need to make that positive
	Rdist=np.clip(Rdist,0,1e12)
	for j in range(1,2):
			Rdist[:,:,:,0,j]=np.min(Rdist[:,:,:,0,:j+1],axis=3)
	for i in range(1,3):
		Rdist[:,:,:,i,0]=np.min(Rdist[:,:,:,i-1:i+1,0],axis=3)
		for j in range(1,2):
			Rdist[:,:,:,i,j]=np.min(np.min(Rdist[:,:,:,:i+1,:j+1],axis=4),axis=3)
	for i in range(3):
		for j in range(2):
			RSum[:,i,j,0]=np.sum(np.sum(Rdist[:,:,:,i,j],axis=1),axis=0)
	RdistN=np.zeros((N0,N1,50,3,2))
	#decreases in #spikes
	for i in range(3):
		RdistN[:,:,:,i,0]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0]=np.max(np.concatenate(\
					(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:,i][:,:,:,None]),axis=3),axis=3)
		RdistN[:,:,:,i,1]=DN[:,:,:,i]
		for k in range(i+1):
			for j in range(i+1):
				RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1]=np.max(np.concatenate(\
				(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2,i][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:,i][:,:,:,None]),axis=3),axis=3)
		#RdistN[:,:,:,i,2]=DN[:,:,:,i]
		#for k in range(i+1):
		#	for j in range(i+1):
		#		RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,2:N2-2,i,2]=np.max(np.concatenate(\
		#		(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,2:N2-2,i,2][:,:,:,None]\
		#		,DN[k:N0-i+k,j:N1-i+j,:-4,i][:,:,:,None]\
		#		,DN[k:N0-i+k,j:N1-i+j,3:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:-2,i][:,:,:,None]\
		#		,DN[k:N0-i+k,j:N1-i+j,1:-3,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,4:,i][:,:,:,None]),axis=3),axis=3)
	#need to make that negative
	RdistN=np.clip(RdistN,-1e12,0)
	for j in range(1,2):
			RdistN[:,:,:,0,j]=np.max(RdistN[:,:,:,0,:j+1],axis=3)
	for i in range(1,3):
		RdistN[:,:,:,i,0]=np.max(RdistN[:,:,:,i-1:i+1,0],axis=3)
		for j in range(1,2):
			RdistN[:,:,:,i,j]=np.max(np.max(RdistN[:,:,:,:i+1,:j+1],axis=4),axis=3)
	for i in range(3):
		for j in range(2):
			RSum[:,i,j,1]=np.sum(np.sum(RdistN[:,:,:,i,j],axis=1),axis=0)
	dD=Rdist[:,:,:,2,1]
	dD1=RdistN[:,:,:,2,1]
	Rgbmap=np.zeros((50,3))
	Rgbmap[5:45,1]=np.bartlett(40)
	Rgbmap[5:25,0]=np.bartlett(40)[20:]
	Rgbmap[:5,0]=np.bartlett(40)[15:20]
	Rgbmap[-5:,0]=np.bartlett(40)[:5]
	Rgbmap[25:-5,2]=np.bartlett(40)[:20]
	Rgbmap[-5:,2]=np.bartlett(40)[20:25]
	Rgbmap[:5,2]=np.bartlett(40)[35:]
	#find difference maxima
	#D=np.cumsum(np.concatenate((Rdist[:,:,:8,0,0][:,:,::-1]\
	#,Rdist[:,:,:,0,0],Rdist[:,:,-7:,0,0][:,:,::-1]),axis=2),axis=2)
	#normalization???-- that one seems reasonable...
	#need to convert that into a complex vector
	A=np.exp(2j*np.pi*np.dot(np.arange(50)/50.,(H1+H2)[:,:,:,None]))[:,:,0]
	AN=np.clip(np.sqrt(np.sum(H1+H2,axis=2)),1,1e12)
	A0=np.cumsum(np.concatenate(((H2+H1)[:,:,:3][:,:,::-1],(H2+H1),(H2+H1)[:,:,-2:][:,:,::-1]),axis=2),axis=2)
	#dD=A[:,:,5:]-A[:,:,:-5]
	DN=np.zeros((N0,N1,3),dtype='complex128')
	DN[:,:,0]=A.copy()
	for i in range(2,4):
		X=np.cumsum(np.concatenate((np.zeros((1,A.shape[1])),A),axis=0),axis=0)
		Y=np.cumsum(np.concatenate((np.zeros((A.shape[0]-i+1,1)),X[i:,:]-X[:-i,:]),axis=1),axis=1)
		DN[(i-1)/2:N0-i/2,(i-1)/2:N1-i/2,i-1]=(Y[:,i:]-Y[:,:-i])*1./i
	RSum2=np.zeros((3),dtype='complex128')
	Rdist2=np.zeros((N0,N1,3),dtype='complex128')
	#increases in #spikes
	for i in range(3):
		Rdist2[:,:,i]=DN[:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					Rdist2[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,i]=np.min(np.abs(np.concatenate(\
					(Rdist2[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,i][:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,i][:,:,None]),axis=2)),axis=2)
	for j in range(1,3):
			Rdist2[:,:,j]=np.min(Rdist2[:,:,:j+1],axis=2)
	dA=np.clip(np.sqrt(A0[:,:,5:]-A0[:,:,:-5]),1,1e12)
	#spatial histo (ranked) of abs and angle
	RdistAbs=np.reshape(np.abs(Rdist2[:,:,2]*1./AN**2).flatten(),(N0,N1))
	RdistAngle=np.reshape(0.5*(np.angle\
	(Rdist2[:,:,2])/np.pi+2*(np.angle(Rdist2[:,:,2])<0)).flatten(),(N0,N1))
	print RdistAbs.max()
	RdistAbsH=np.array(50*np.argsort(np.argsort(RdistAbs.flatten()))*1./(N0*N1),dtype=int)
	RdistAngleH=np.array(49.999*RdistAngle,dtype=int)
	#difference map
	W=np.clip(np.dot(np.reshape((dD/dA),(-1,50)),Rgbmap-0.5)/50.+0.5,0,1)
	W1=np.clip(np.dot(np.reshape((dD1/dA),(-1,50)),Rgbmap-0.5)/50.+0.5,0,1)
	#sparse maps (want abs and angle --> extra dimension)
	WX=np.zeros((N0*N1,50))
	WY=np.zeros((N0*N1,50))
	WNX=np.zeros((N0*N1,50))
	WNY=np.zeros((N0*N1,50))
	WX[np.arange(N0*N1,dtype=int),RdistAbsH.flatten()]=(AN.flatten())#RdistAngleH.flatten()
	WY[np.arange(N0*N1,dtype=int),RdistAngleH.flatten()]=RdistAbsH.flatten()
	WNX[np.arange(N0*N1,dtype=int),RdistAbsH.flatten()]=(AN.flatten())
	WNY[np.arange(N0*N1,dtype=int),RdistAngleH.flatten()]=1.#(AN.flatten())
	#sum to obtain histogram
	WRX=np.einsum('ij,ik',WX[:,:],np.clip((np.reshape((dD/dA),(-1,50)))/50.,0.,1.))
	WRY=np.einsum('ij,ik',WY[:,:],np.clip((np.reshape((dD/dA),(-1,50)))/50.,0.,1.))
	WRNX=np.einsum('ij,ik',WNX[:,:],np.clip((np.reshape((dD/dA),(-1,50)))/50.,0.,1.))
	WRNY=np.einsum('ij,ik',WNY[:,:],np.clip((np.reshape((dD/dA),(-1,50)))/50.,0.,1.))
	WRXl=np.einsum('ij,ik',WX[:,:],np.clip((np.reshape((-dD1/dA),(-1,50)))/50.,0.,1.))
	WRYl=np.einsum('ij,ik',WY[:,:],np.clip((np.reshape((-dD1/dA),(-1,50)))/50.,0.,1.))
	WRNXl=np.einsum('ij,ik',WNX[:,:],np.clip((np.reshape((-dD1/dA),(-1,50)))/50.,0.,1.))
	WRNYl=np.einsum('ij,ik',WNY[:,:],np.clip((np.reshape((-dD1/dA),(-1,50)))/50.,0.,1.))
	'''
	WRY=np.sum(WY[:,None,:]*np.reshape(Rdist[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	WRNX=np.sum(WNX[:,None,:]*np.reshape(Rdist[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	WRNY=np.sum(WNY[:,None,:]*np.reshape(Rdist[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	WRXl=np.sum(WX[:,None,:]*np.reshape(RdistN[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	WRYl=np.sum(WY[:,None,:]*np.reshape(RdistN[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	WRNXl=np.sum(WNX[:,None,:]*np.reshape(RdistN[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	WRNYl=np.sum(WNY[:,None,:]*np.reshape(RdistN[:,:,:,2,1],(-1,50))[:,:,None],axis=0)
	'''
	#convert into colour
	W2=np.clip(np.sqrt(WRNX.flatten())[:,None]\
	*(Rgbmap[np.array((WRX*1./np.clip(WRNX,1,1e12)).flatten(),dtype=int),:]-0.33)/10.+0.33,0,1)
	W3=np.clip(np.sqrt(WRNY.flatten())[:,None]\
	*(Rgbmap[np.array((WRY*1./np.clip(WRNY,1,1e12)).flatten(),dtype=int),:]-0.33)/10.+0.33,0,1)
	W2l=np.clip(np.sqrt(WRNXl.flatten())[:,None]\
	*(Rgbmap[np.array((WRXl*1./np.clip(WRNXl,1,1e12)).flatten(),dtype=int),:]-0.33)/10.+0.33,0,1)
	W3l=np.clip(np.sqrt(WRNYl.flatten())[:,None]\
	*(Rgbmap[np.array((WRYl*1./np.clip(WRNYl,1,1e12)).flatten(),dtype=int),:]-0.33)/10.+0.33,0,1)
	#print W2.max(), W2.min(), W3.max(), W3.min(), W2l.max(), W2l.min(), W3l.max(), W3l.min()
	W2=np.reshape(W2,(50,50,3))
	W2l=np.reshape(W2l,(50,50,3))
	W3=np.reshape(W3,(50,50,3))
	W3l=np.reshape(W3l,(50,50,3))
	W23=np.concatenate((np.concatenate((W2,W3),axis=1),np.concatenate((W2l,W3l),axis=1)),axis=0)
	#for plotting angle and abs#
	W4=np.clip(np.dot(WX*1./np.clip(WNX,1,1e12),Rgbmap-0.5)*np.sqrt(np.sum(WNX,axis=1))[:,None]/10.+0.5,0,1)
	W5=np.clip(np.dot(WY*1./np.clip(WNY,1,1e12),Rgbmap-0.5)*np.sqrt(np.sum(WNY,axis=1))[:,None]/10.+0.5,0,1)
	#W3=np.clip(np.sqrt(WRN.flatten())[:,None]*(Rgbmap[RdistAngleH.flatten()]-0.33)/10.+0.33,0,1)
	#WR=np.einsum('i,jk',RdistAbsH,W)

	
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (1./10., 1., 1.), (5./10., 0.0, 0.0),\
	(9./10., 0.0, 0.0), (1.0, 0.25, 0.25)),\
	'green': ((0.0,0.,0.),(1./10.,0.,0.),(5./10.,1.,1.),(9./10.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(1./10.,0.,0.),(5./10.,0.,0.),(9./10.,1.,1.),(1.0, 0.75, 0.75))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (1./10., 0.0, 0.0), (5./10., 1., 1.),\
	(9./10., 1., 1.), (1.0, 0.75, 0.75)),\
	'green': ((0.0,1.,1.),(1./10.,1.,1.),(5./10.,0.0,0.0),(9./10.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(1./10.,1.,1.),(5./10.,1.,1.),(9./10.,0.0,0.0),(1.0, 0.25, 0.25))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (low)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	#cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (high)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	#ax.imshow((Q*1./(QN+(QN==0))[:,None]).transpose()\
	#ax.imshow(S.transpose()\
	#,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	#,interpolation='none',origin='lower',extent=(-100,100,0,100))
	#ax.set_ylim(0,100)
	#ax.set_xlim(-100,100)
	#ax.set_xticks(np.array([-100,0,100]))
	#ax.set_yticks(np.array([0,50,100]))np.reshape(W,((yMax-yMin)*Res,-1,3))
	ax.imshow(np.reshape(W4,((yMax-yMin)*Res,-1,3)), aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(np.reshape(W5,((yMax-yMin)*Res,-1,3))\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	'''
	pylab.figure()
	Clr=np.array(['k','b','r'])
	Clr2=np.array(['y','c','m'])
	L=np.array(['-','--',':'])
	for i in range(3):
		for j in range(2):
			pylab.plot(np.arange(50),RSum[:,i,j,0],Clr[i]+L[j])
			pylab.plot(np.arange(50),-RSum[:,i,j,1],Clr2[i]+L[j])
	'''
	#c=np.histogram2d(np.sign(b)*50+a,np.sign(b1)*50+a1,bins=(np.arange(-50,150),np.arange(-50,150)))[0]
	fig2=pylab.figure(1,figsize=(8,8))
	#fig2.text(0.05,0.5,'typical excitability rank',rotation='vertical',ha='center',va='center')#,fontsize='xx-large'
	#fig2.text(0.5,0.05,'typical excitability change',ha='center',va='center')#,fontsize='xx-large'
	#fig2.text(0.97,0.3,'neg.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	#fig2.text(0.97,0.7,'pos.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	#fig2.text(0.3,0.97,'neg. change',ha='center',va='center',fontsize='small')
	#fig2.text(0.7,0.97,'pos. change',ha='center',va='center',fontsize='small')
	ax4=fig2.add_axes([0.1, 0.1, 0.8, 0.8],axisbg=bgColor)
	ax5=fig2.add_axes([0.1, 0.9, 0.4, 0.03])
	ax6=fig2.add_axes([0.5, 0.9, 0.4, 0.03])
	ax7=fig2.add_axes([0.9, 0.1, 0.03, 0.4])
	ax8=fig2.add_axes([0.9, 0.5, 0.03, 0.4])
	ax4.set_xticks(np.array([]))
	ax4.set_yticks(np.array([]))
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	ax6.set_xticks(np.array([]))
	ax6.set_yticks(np.array([]))
	ax7.set_xticks(np.array([]))
	ax7.set_yticks(np.array([]))
	ax8.set_xticks(np.array([]))
	ax8.set_yticks(np.array([]))
	cbar5 =mplt.colorbar.ColorbarBase(ax=ax5,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar6 =mplt.colorbar.ColorbarBase(ax=ax6,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar7 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar8 =mplt.colorbar.ColorbarBase(ax=ax8,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar5.set_ticks(np.array([]))
	cbar6.set_ticks(np.array([]))
	cbar7.set_ticks(np.array([]))
	cbar8.set_ticks(np.array([]))
	ax4.imshow(W23, aspect='equal'\
	, interpolation='none', origin='lower', extent=(0,200,0,200))
	ax4.set_ylim(0,200)
	ax4.set_xlim(0,200)
	#ax4.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	#ax4.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	ax4.plot(np.array([100,100]),np.array([0,200]),'k-')
	ax4.plot(np.array([0,200]),np.array([100,100]),'k-')
	pylab.savefig(FileName2)
	pylab.close()
	#pylab.figure()
	#pylab.imshow(np.reshape(Hsum,((yMax-yMin)*Res,-1))\
	#,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	#pylab.colorbar()
	#pylab.show()
	return

#construct a map for ranking the sum of excitability at each location.
def TransitionMatricesOpt3(HdfFile1, HdfFile2, FileName, FileName2\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(50.*np.exp(-cQ1)))*50./len(cQ1)
	A2=np.argsort(np.argsort(50.*np.exp(-cQ2)))*50./len(cQ2)
	L1=len(A1)
	L2=len(A2)
	A3=np.zeros((L1,4))
	A4=np.zeros((L2,4))
	for i in range(4):
		Az=np.argsort(np.argsort(np.concatenate((A1+1e-6*scipy.rand(L1)+1e3*np.array((Res/2)\
		*((Loc1[:,0]-(i%2)*1./Res)+64*Res/2*(Loc1[:,1]-(i/2)*1./Res)),dtype=int)\
		,A2+1e-6*scipy.rand(L2)+1e3*np.array((Res/2)\
		*((Loc2[:,0]-(i%2)*1./Res)+64*Res/2*(Loc2[:,1]-(i/2)*1./Res)),dtype=int))))) #--> not that easy
		Aind=np.concatenate(\
		(Res/2*(np.array(np.clip(Loc1[:,0]-(i%2)*1./Res,0,1e6),dtype=int)\
		+64*Res/2*(np.array(np.clip(Loc1[:,1]-(i/2)*1./Res,0,1e6),dtype=int)))\
		,Res/2*(np.array(np.clip(Loc2[:,0]-(i%2)*1./Res,0,1e6),dtype=int)\
		+64*Res/2*(np.array(np.clip(Loc2[:,1]-(i/2)*1./Res,0,1e6),dtype=int)))))
		Ah=np.concatenate((np.array([0]),np.cumsum(np.histogram(Aind\
		,bins=np.arange((Res/2)**2*4096+1))[0])))
		print Ah.shape, Az.shape, np.sum(np.diff(Ah)), Ah[-1], Loc1.min(), Loc2.min()
		Ar=np.argsort(Az)
		Az[Ar]-=np.repeat(Ah[:-1],np.diff(Ah))#Ah[Aind]#
		Az[Ar]*=50./np.clip(np.repeat(np.diff(Ah),np.diff(Ah)),1,1e12)#(1./np.clip(np.diff(Ah),1,1e12))[Aind]
		A3[:,i]=Az[:L1]
		A4[:,i]=Az[L1:]
	print A3.max(),A3.min()
	#could repeat that procedure with a shifted version and average/minimize
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(51)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(51)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	D=np.cumsum(np.concatenate(((H2-H1)[:,:,:3][:,:,::-1],(H2-H1),(H2-H1)[:,:,-2:][:,:,::-1]),axis=2),axis=2)
	dD=D[:,:,5:]-D[:,:,:-5]
	N0=dD.shape[0]
	N1=dD.shape[1]
	N2=dD.shape[2]
	#-----------------
	H3=np.zeros((N0+1,N1+1,N2,4))
	H4=np.zeros((N0+1,N1+1,N2,4))
	for i in range(4):
		H3[:,:,:,i]=np.histogramdd(np.concatenate((Loc1,A3[:,i][:,None]),axis=1)\
		,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
		,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(51)))[0]
		H4[:,:,:,i]=np.histogramdd(np.concatenate((Loc2,A4[:,i][:,None]),axis=1)\
		,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
		,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(51)))[0]
	H3=H3[1:,1:,:,:]+H3[:-1,1:,:,:]+H3[:-1,:-1,:,:]+H3[1:,:-1,:,:]
	H4=H4[1:,1:,:,:]+H4[:-1,1:,:,:]+H4[:-1,:-1,:,:]+H4[1:,:-1,:,:]
	D3=np.cumsum(np.concatenate(((H4-H3)[:,:,:3,:][:,:,::-1,:],(H4-H3),(H4-H3)[:,:,-2:,:][:,:,::-1,:]),axis=2),axis=2)
	dD3=D3[:,:,5:,:]-D3[:,:,:-5,:]
	#min absolute value
	D3min=np.clip(np.min(dD3,axis=3),0,1e12)+np.clip(np.max(dD3,axis=3),-1e12,0)
	#-----------------------
	DN=np.zeros((N0,N1,N2,3))
	DN[:,:,:,0]=dD.copy()
	for i in range(2,4):
		X=np.cumsum(np.concatenate((np.zeros((1,dD.shape[1],dD.shape[2])),dD),axis=0),axis=0)
		Y=np.cumsum(np.concatenate((np.zeros((dD.shape[0]-i+1,1,dD.shape[2])),X[i:,:,:]-X[:-i,:,:]),axis=1),axis=1)
		DN[(i-1)/2:N0-i/2,(i-1)/2:N1-i/2,:,i-1]=(Y[:,i:,:]-Y[:,:-i,:])*1./i
	RSum=np.zeros((50,3,2,2))
	Rdist=np.zeros((N0,N1,50,3,2))
	#increases in #spikes
	for i in range(3):
		Rdist[:,:,:,i,0]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0]=np.min(np.concatenate(\
					(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:,i][:,:,:,None]),axis=3),axis=3)
		Rdist[:,:,:,i,1]=DN[:,:,:,i]
		for k in range(i+1):
			for j in range(i+1):
				Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1]=np.min(np.concatenate(\
				(Rdist[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2,i][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:,i][:,:,:,None]),axis=3),axis=3)
	#need to make that positive
	Rdist=np.clip(Rdist,0,1e12)
	for j in range(1,2):
		Rdist[:,:,:,0,j]=np.min(Rdist[:,:,:,0,:j+1],axis=3)
	for i in range(1,3):
		Rdist[:,:,:,i,0]=np.min(Rdist[:,:,:,i-1:i+1,0],axis=3)
		for j in range(1,2):
			Rdist[:,:,:,i,j]=np.min(np.min(Rdist[:,:,:,:i+1,:j+1],axis=4),axis=3)
	for i in range(3):
		for j in range(2):
			RSum[:,i,j,0]=np.sum(np.sum(Rdist[:,:,:,i,j],axis=1),axis=0)
	RdistN=np.zeros((N0,N1,50,3,2))
	#decreases in #spikes
	for i in range(3):
		RdistN[:,:,:,i,0]=DN[:,:,:,i]
		if i:
			for k in range(i+1):
				for j in range(i+1):
					RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0]=np.max(np.concatenate(\
					(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,:,i,0][:,:,:,None]\
					,DN[k:N0-i+k,j:N1-i+j,:,i][:,:,:,None]),axis=3),axis=3)
		RdistN[:,:,:,i,1]=DN[:,:,:,i]
		for k in range(i+1):
			for j in range(i+1):
				RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1]=np.max(np.concatenate(\
				(RdistN[i/2:N0-(i+1)/2,i/2:N1-(i+1)/2,1:N2-1,i,1][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,:-2,i][:,:,:,None]\
				,DN[k:N0-i+k,j:N1-i+j,1:-1,i][:,:,:,None],DN[k:N0-i+k,j:N1-i+j,2:,i][:,:,:,None]),axis=3),axis=3)
	#need to make that negative
	RdistN=np.clip(RdistN,-1e12,0)
	for j in range(1,2):
			RdistN[:,:,:,0,j]=np.max(RdistN[:,:,:,0,:j+1],axis=3)
	for i in range(1,3):
		RdistN[:,:,:,i,0]=np.max(RdistN[:,:,:,i-1:i+1,0],axis=3)
		for j in range(1,2):
			RdistN[:,:,:,i,j]=np.max(np.max(RdistN[:,:,:,:i+1,:j+1],axis=4),axis=3)
	for i in range(3):
		for j in range(2):
			RSum[:,i,j,1]=np.sum(np.sum(RdistN[:,:,:,i,j],axis=1),axis=0)
	dD=Rdist[:,:,:,2,1]
	dD1=RdistN[:,:,:,2,1]
	Rgbmap=np.zeros((50,3))
	Rgbmap[5:45,1]=np.bartlett(40)
	Rgbmap[5:25,0]=np.bartlett(40)[20:]
	Rgbmap[:5,0]=np.bartlett(40)[15:20]
	Rgbmap[-5:,0]=np.bartlett(40)[:5]
	Rgbmap[25:-5,2]=np.bartlett(40)[:20]
	Rgbmap[-5:,2]=np.bartlett(40)[20:25]
	Rgbmap[:5,2]=np.bartlett(40)[35:]
	#make a histogram using the outer product and subtract small changes
	# by removing the outer product of the histogram intersection (smoothened HI eventually)
	A=np.cumsum(np.concatenate(((H2+H1)[:,:,:3][:,:,::-1],(H2+H1),(H2+H1)[:,:,-2:][:,:,::-1]),axis=2),axis=2)
	N=np.clip(A[:,:,-1],0,1e3)#clipping might be useful for normalization (maximum weight)
	dA=np.clip(np.sqrt(A[:,:,5:]-A[:,:,:-5]),1,1e12)
	Anp=np.zeros((100,100))
	#HInp=np.zeros((100,100))
	for i in range(N0*N1):
		xa=np.concatenate((-dD1[i%N0,i/N0,:],dD[i%N0,i/N0,:]))#all positive values
		xa*=1./np.clip(np.sum(xa),1,1e12)
		x0=np.concatenate((-np.clip(D3min[i%N0,i/N0,:],-1e12,0),np.clip(D3min[i%N0,i/N0,:],0,1e12)))
		x0*=1./np.clip(np.sum(x0),1,1e12)
		HI=np.min(np.concatenate((xa[:,None],x0[:,None]),axis=1),axis=1)
		Anp+=N[i%N0,i/N0]*(np.outer(xa,x0))#-np.outer(HI,HI))#weighting? -- number of spikes (due to ranking)
	#difference map
	W=np.clip(np.dot(np.reshape((dD/dA),(-1,50)),Rgbmap-0.5)/50.+0.5,0,1)
	W1=np.clip(np.dot(np.reshape((dD1/dA),(-1,50)),Rgbmap-0.5)/50.+0.5,0,1)
	W2=np.clip(np.dot(np.reshape((D3min/dA),(-1,50)),Rgbmap-0.5)/200.+0.5,0,1)
	W3=np.clip((W1+W-W2),0,1)
	
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (1./10., 1., 1.), (5./10., 0.0, 0.0),\
	(9./10., 0.0, 0.0), (1.0, 0.25, 0.25)),\
	'green': ((0.0,0.,0.),(1./10.,0.,0.),(5./10.,1.,1.),(9./10.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(1./10.,0.,0.),(5./10.,0.,0.),(9./10.,1.,1.),(1.0, 0.75, 0.75))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (1./10., 0.0, 0.0), (5./10., 1., 1.),\
	(9./10., 1., 1.), (1.0, 0.75, 0.75)),\
	'green': ((0.0,1.,1.),(1./10.,1.,1.),(5./10.,0.0,0.0),(9./10.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(1./10.,1.,1.),(5./10.,1.,1.),(9./10.,0.0,0.0),(1.0, 0.25, 0.25))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (low)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	#cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (high)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	#ax.imshow((Q*1./(QN+(QN==0))[:,None]).transpose()\
	#ax.imshow(S.transpose()\
	#,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	#,interpolation='none',origin='lower',extent=(-100,100,0,100))
	#ax.set_ylim(0,100)
	#ax.set_xlim(-100,100)
	#ax.set_xticks(np.array([-100,0,100]))
	#ax.set_yticks(np.array([0,50,100]))np.reshape(W,((yMax-yMin)*Res,-1,3))
	ax.imshow(np.reshape((W+W1)/2.,((yMax-yMin)*Res,-1,3)), aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(np.reshape(W3,((yMax-yMin)*Res,-1,3))\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	'''
	pylab.figure()
	Clr=np.array(['k','b','r'])
	Clr2=np.array(['y','c','m'])
	L=np.array(['-','--',':'])
	for i in range(3):
		for j in range(2):
			pylab.plot(np.arange(50),RSum[:,i,j,0],Clr[i]+L[j])
			pylab.plot(np.arange(50),-RSum[:,i,j,1],Clr2[i]+L[j])
	'''
	#c=np.histogram2d(np.sign(b)*50+a,np.sign(b1)*50+a1,bins=(np.arange(-50,150),np.arange(-50,150)))[0]
	fig2=pylab.figure(1,figsize=(8,8))
	#fig2.text(0.05,0.5,'typical excitability rank',rotation='vertical',ha='center',va='center')#,fontsize='xx-large'
	#fig2.text(0.5,0.05,'typical excitability change',ha='center',va='center')#,fontsize='xx-large'
	#fig2.text(0.97,0.3,'neg.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	#fig2.text(0.97,0.7,'pos.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	#fig2.text(0.3,0.97,'neg. change',ha='center',va='center',fontsize='small')
	#fig2.text(0.7,0.97,'pos. change',ha='center',va='center',fontsize='small')
	ax4=fig2.add_axes([0.1, 0.1, 0.8, 0.8],axisbg=bgColor)
	ax5=fig2.add_axes([0.1, 0.9, 0.4, 0.03])
	ax6=fig2.add_axes([0.5, 0.9, 0.4, 0.03])
	ax7=fig2.add_axes([0.9, 0.1, 0.03, 0.4])
	ax8=fig2.add_axes([0.9, 0.5, 0.03, 0.4])
	ax4.set_xticks(np.array([]))
	ax4.set_yticks(np.array([]))
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	ax6.set_xticks(np.array([]))
	ax6.set_yticks(np.array([]))
	ax7.set_xticks(np.array([]))
	ax7.set_yticks(np.array([]))
	ax8.set_xticks(np.array([]))
	ax8.set_yticks(np.array([]))
	cbar5 =mplt.colorbar.ColorbarBase(ax=ax5,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar6 =mplt.colorbar.ColorbarBase(ax=ax6,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar7 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar8 =mplt.colorbar.ColorbarBase(ax=ax8,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar5.set_ticks(np.array([]))
	cbar6.set_ticks(np.array([]))
	cbar7.set_ticks(np.array([]))
	cbar8.set_ticks(np.array([]))
	ax4.imshow(np.log10(np.clip(Anp,1,1e4)), aspect='equal', vmin=0., vmax=4., cmap=SpikesCmap\
	, interpolation='none', origin='lower', extent=(0,200,0,200))
	ax4.set_ylim(0,200)
	ax4.set_xlim(0,200)
	#ax4.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	#ax4.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	ax4.plot(np.array([100,100]),np.array([0,200]),'k-')
	ax4.plot(np.array([0,200]),np.array([100,100]),'k-')
	pylab.savefig(FileName2)
	pylab.close()
	#pylab.figure()
	#pylab.imshow(np.reshape(Hsum,((yMax-yMin)*Res,-1))\
	#,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	#pylab.colorbar()
	#pylab.show()
	return

def TransitionMatricesPlotDev2(HdfFile1, HdfFile2, HdfFile3, FileName, FileName2\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	g=h5py.File(HdfFile3,'r')
	Ind3=g[IndFieldName].value
	Amp3=g[AmpFieldName].value
	Amp3=Amp3[Ind3]
	cQ3=g[FieldName].value
	cQ3=cQ3[Ind3]
	pAmpInd3=(Amp3>=ampThreshold)*(cQ3>=cThreshold)
	cQ3=cQ3[pAmpInd3]
	Loc3=g[LocFieldName].value[Ind3,:][pAmpInd3,:]
	Amp3=Amp3[pAmpInd3]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(100.*np.exp(-cQ1)))*100./len(cQ1)
	A2=np.argsort(np.argsort(100.*np.exp(-cQ2)))*100./len(cQ2)
	A3=np.argsort(np.argsort(100.*np.exp(-cQ3)))*100./len(cQ3)
	L1=len(A1)
	L2=len(A2)
	L3=len(A3)
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H3=np.histogramdd(np.concatenate((Loc3,A3[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	H3=H3[1:,1:,:]+H3[:-1,1:,:]+H3[:-1,:-1,:]+H3[1:,:-1,:]
	Rgbmap=np.zeros((100,3))
	Rgbmap[10:90,1]=np.bartlett(80)
	Rgbmap[10:50,0]=np.bartlett(80)[40:]
	Rgbmap[:10,0]=np.bartlett(80)[30:40]
	Rgbmap[-10:,0]=np.bartlett(80)[:10]
	Rgbmap[50:-10,2]=np.bartlett(80)[:40]
	Rgbmap[-10:,2]=np.bartlett(80)[40:50]
	Rgbmap[:10,2]=np.bartlett(80)[70:]
	#find difference maxima
	D=np.cumsum(np.concatenate(((H2-H1)[:,:,:6][:,:,::-1],(H2-H1),(H2-H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	D1=np.cumsum(np.concatenate(((H3-H2)[:,:,:6][:,:,::-1],(H3-H2),(H3-H2)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#need to subtract a uniform distribution here...
	A=np.cumsum(np.concatenate(((H2+H1)[:,:,:6][:,:,::-1],(H2+H1),(H2+H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	A1=np.cumsum(np.concatenate(((H3+H2)[:,:,:6][:,:,::-1],(H3+H2),(H3+H2)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#print D.shape, H1.shape, H2.shape
	dD=D[:,:,11:]-D[:,:,:-11]
	dD1=D1[:,:,11:]-D1[:,:,:-11]
	dA=np.clip(np.sqrt(A[:,:,11:]-A[:,:,:-11]),1,1e12)
	dA1=np.clip(np.sqrt(A1[:,:,11:]-A1[:,:,:-11]),1,1e12)
	a=np.argmax(np.abs(dD/dA),axis=2).flatten()
	a1=np.argmax(np.abs(dD1/dA1),axis=2).flatten()
	b=np.clip(((dD/dA).flatten()[a+np.arange(dD.shape[0]*dD.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	b1=np.clip(((dD1/dA1).flatten()[a1+np.arange(dD1.shape[0]*dD1.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	#convert to color
	X=np.reshape(Rgbmap[a,:],((yMax-yMin)*Res,-1,3))
	X1=np.reshape(Rgbmap[a1,:],((yMax-yMin)*Res,-1,3))
	Y=(X-0.5)*np.reshape(2*b,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	Y1=(X1-0.5)*np.reshape(2*b1,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	
	#W=np.clip(np.dot(np.reshape((dD/dA),(-1,100))-0.5,Rgbmap-0.5)/150.+0.5,0,1)
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (1./10., 1., 1.), (5./10., 0.0, 0.0),\
	(9./10., 0.0, 0.0), (1.0, 0.25, 0.25)),\
	'green': ((0.0,0.,0.),(1./10.,0.,0.),(5./10.,1.,1.),(9./10.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(1./10.,0.,0.),(5./10.,0.,0.),(9./10.,1.,1.),(1.0, 0.75, 0.75))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (1./10., 0.0, 0.0), (5./10., 1., 1.),\
	(9./10., 1., 1.), (1.0, 0.75, 0.75)),\
	'green': ((0.0,1.,1.),(1./10.,1.,1.),(5./10.,0.0,0.0),(9./10.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(1./10.,1.,1.),(5./10.,1.,1.),(9./10.,0.0,0.0),(1.0, 0.25, 0.25))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	
	
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (decreasing)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	#cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (increasing)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	ax.imshow(Y1, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(Y\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	c=np.histogram2d(np.sign(b)*50+a,np.sign(b1)*50+a1,bins=(np.arange(-50,150),np.arange(-50,150)))[0]
	fig2=pylab.figure(1,figsize=(8,8))
	fig2.text(0.05,0.5,'typical excitability change pre',rotation='vertical',ha='center',va='center')#,fontsize='xx-large'
	fig2.text(0.5,0.05,'typical excitability change post',ha='center',va='center')#,fontsize='xx-large'
	fig2.text(0.97,0.3,'neg. change',rotation='vertical',fontsize='small',ha='center',va='center')
	fig2.text(0.97,0.7,'pos. change',rotation='vertical',fontsize='small',ha='center',va='center')
	fig2.text(0.3,0.97,'neg. change',ha='center',va='center',fontsize='small')
	fig2.text(0.7,0.97,'pos. change',ha='center',va='center',fontsize='small')
	ax4=fig2.add_axes([0.1, 0.1, 0.8, 0.8],axisbg=bgColor)
	ax5=fig2.add_axes([0.1, 0.9, 0.4, 0.03])
	ax6=fig2.add_axes([0.5, 0.9, 0.4, 0.03])
	ax7=fig2.add_axes([0.9, 0.1, 0.03, 0.4])
	ax8=fig2.add_axes([0.9, 0.5, 0.03, 0.4])
	ax4.set_xticks(np.array([]))
	ax4.set_yticks(np.array([]))
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	ax6.set_xticks(np.array([]))
	ax6.set_yticks(np.array([]))
	ax7.set_xticks(np.array([]))
	ax7.set_yticks(np.array([]))
	ax8.set_xticks(np.array([]))
	ax8.set_yticks(np.array([]))
	cbar5 =mplt.colorbar.ColorbarBase(ax=ax5,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar6 =mplt.colorbar.ColorbarBase(ax=ax6,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar7 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar8 =mplt.colorbar.ColorbarBase(ax=ax8,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar5.set_ticks(np.array([]))
	cbar6.set_ticks(np.array([]))
	cbar7.set_ticks(np.array([]))
	cbar8.set_ticks(np.array([]))
	ax4.imshow(np.log10(np.clip(c.transpose(),1,1e3)), aspect='equal', vmin=0., vmax=3., cmap=SpikesCmap\
	, interpolation='none', origin='lower', extent=(0,200,0,200))
	ax4.set_ylim(0,200)
	ax4.set_xlim(0,200)
	ax4.plot(np.array([100,100]),np.array([0,200]),'k-')
	ax4.plot(np.array([0,200]),np.array([100,100]),'k-')
	pylab.savefig(FileName2)
	pylab.close()
	return


def TransitionMatricesPlotFiringrate(HdfFile1, HdfFile2, FileName, FileName2\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(15,7))
	ax=fig.add_axes([0.04, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax1=fig.add_axes([0.81, 0.13, 0.015, 0.8])
	ax2=fig.add_axes([0.44, 0.13, 0.36, 0.8],axisbg=bgColor)
	ax3=fig.add_axes([0.89, 0.13, 0.015, 0.8])
	ax.tick_params(pad=8)
	ax2.tick_params(pad=8)
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	#Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	#Amp1=Amp1[Ind1]
	#cQ1=g[FieldName].value
	#cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)#*(cQ1>=cThreshold)
	#cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	#Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	#Amp2=Amp2[Ind2]
	#cQ2=g[FieldName].value
	#cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)#*(cQ2>=cThreshold)
	#cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	#histogram of spike count
	HN1=np.histogram2d(Loc1[:,0],Loc1[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN2=np.histogram2d(Loc2[:,0],Loc2[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN1=HN1[1:,1:]+HN1[:-1,1:]+HN1[:-1,:-1]+HN1[1:,:-1]
	HN2=HN2[1:,1:]+HN2[:-1,1:]+HN2[:-1,:-1]+HN2[1:,:-1]
	#cQ1=HN1[np.array((Loc1[:,0]+yMin)*Res,dtype=int),np.array((Loc1[:,1]+xMin)*Res,dtype=int)]
	#cQ1=HN1[np.array((Loc2[:,0]+yMin)*Res,dtype=int),np.array((Loc2[:,1]+xMin)*Res,dtype=int)]
	HN1=HN1.flatten()
	HN2=HN2.flatten()
	A1=np.argsort(np.argsort(100.*np.exp(-HN1)))*100./len(HN1)
	A2=np.argsort(np.argsort(100.*np.exp(-HN2)))*100./len(HN2)
	A3=np.argsort(np.argsort(A1-A2))*100./len(HN2)
	#L1=len(A1)
	#L2=len(A2)
	#histogram of excitability of different locations
	'''
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	'''
	#H1r=np.reshape(H1,(-1,100))
	#H2r=np.reshape(H2,(-1,100))

	#threshold spike count
	#NInd=(HN1+HN2)>=100
	#NInd=np.argsort(HN1+HN2)[-500:]
	#Hsum=np.clip(np.sum((H2r-H1r)*np.arange(100)[None,:],axis=1)*1./np.clip((HN2+HN1),1,1e12),-100,100)
	##Hloc=np.histogram(Hsum,bins=np.arange(-100,101))[0]
	#Histogram of average and single location activity change
	#print Hsum.shape
	#Q=np.histogram2d((Hsum[:,None]*np.ones(100)[None,:]).flatten()\
	#,(np.ones(Hsum.shape[0])[:,None]*np.arange(100)[None,:]).flatten()\
	#,bins=(np.arange(-100,101),np.arange(101)),weights=(H2r-H1r).flatten())[0]
	#QN=np.histogram(Hsum,bins=(np.arange(-100,101)),weights=np.sqrt(np.sum(H2r+H1r,axis=1)))[0]#
	#need to figure out where the largest difference is (in Ex.) take 10% window (sliding, periodic boundary)
	#use absolute difference (^1/3?)
	#cannot plot that...clustering? -- distances? - correlations... or firing rate cutoff take highest 500?
	#S=((H2r-H1r)*1./np.clip((np.sqrt(np.sum(H2r+H1r,axis=1))),1,1e12)[:,None])[NInd,:][np.argsort(Hsum[NInd]),:]
	
	
	Rgbmap=np.zeros((100,3))
	Rgbmap[10:90,1]=np.bartlett(80)
	Rgbmap[10:50,0]=np.bartlett(80)[40:]
	Rgbmap[:10,0]=np.bartlett(80)[30:40]
	Rgbmap[-10:,0]=np.bartlett(80)[:10]
	Rgbmap[50:-10,2]=np.bartlett(80)[:40]
	Rgbmap[-10:,2]=np.bartlett(80)[40:50]
	Rgbmap[:10,2]=np.bartlett(80)[70:]
	#find difference maxima
	#D=np.cumsum(np.concatenate(((H2-H1)[:,:,:6][:,:,::-1],(H2-H1),(H2-H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#need to subtract a uniform distribution here...
	#D1=np.cumsum(np.concatenate(((H1-1.5*np.mean(H1,axis=2)[:,:,None])[:,:,:6][:,:,::-1]\
	#,(H1-1.5*np.mean(H1,axis=2)[:,:,None]),(H1-1.5*np.mean(H1,axis=2)[:,:,None])[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#A=np.cumsum(np.concatenate(((H2+H1)[:,:,:6][:,:,::-1],(H2+H1),(H2+H1)[:,:,-5:][:,:,::-1]),axis=2),axis=2)
	#print D.shape, H1.shape, H2.shape
	#dD=D[:,:,11:]-D[:,:,:-11]
	#dD1=D1[:,:,11:]-D1[:,:,:-11]
	#dA=np.clip(np.sqrt(A[:,:,11:]-A[:,:,:-11]),1,1e12)
	#a=np.argmax(np.abs(dD/dA),axis=2).flatten()
	#a1=np.argmax(np.abs(dD1),axis=2).flatten()
	#a1=np.argsort(np.argsort(np.abs(dD1),axis=2),axis=2).flatten()\
	#[a+np.arange(dD.shape[0]*dD.shape[1],dtype=int)*100]
	#b=np.zeros(a.shape)
	#b=np.clip(((dD/dA).flatten()[a+np.arange(dD.shape[0]*dD.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	#b1=np.clip(((np.sign(dD1)*np.sqrt(np.abs(dD1))).flatten()[a1+np.arange(dD1.shape[0]*dD1.shape[1],dtype=int)*100])/6.,-0.5,0.5)
	#convert to color
	a=np.array(np.abs(2*A3-100.+1e-12),dtype=int)
	a1=np.array(np.abs(2*A1-100.+1e-12),dtype=int)
	b=np.sign(2*A3-100.)
	b1=np.sign(A1-50)
	X=np.reshape(Rgbmap[a,:],((yMax-yMin)*Res,-1,3))
	X1=np.reshape(Rgbmap[a1,:],((yMax-yMin)*Res,-1,3))
	Y=(X-0.5)*np.reshape(2*b,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	Y1=(X1-0.5)*np.reshape(2*b1,((yMax-yMin)*Res,-1))[:,:,None]+0.5
	#(H2-H1)/(H2+H1)
	
	#W=np.clip(np.dot(np.reshape((dD/dA),(-1,100))-0.5,Rgbmap-0.5)/150.+0.5,0,1)
	rgbkdict = {'red': ((0.0, 0.75, 0.75), (1./10., 1., 1.), (5./10., 0.0, 0.0),\
	(9./10., 0.0, 0.0), (1.0, 0.25, 0.25)),\
	'green': ((0.0,0.,0.),(1./10.,0.,0.),(5./10.,1.,1.),(9./10.,0.,0.),(1.0, 0., 0.)),\
	'blue': ((0.0,0.25,0.25),(1./10.,0.,0.),(5./10.,0.,0.),(9./10.,1.,1.),(1.0, 0.75, 0.75))}
	rgbwdict = {'red': ((0.0, 0.25, 0.25), (1./10., 0.0, 0.0), (5./10., 1., 1.),\
	(9./10., 1., 1.), (1.0, 0.75, 0.75)),\
	'green': ((0.0,1.,1.),(1./10.,1.,1.),(5./10.,0.0,0.0),(9./10.,1.,1.),(1.0, 1., 1.)),\
	'blue': ((0.0,0.75,0.75),(1./10.,1.,1.),(5./10.,1.,1.),(9./10.,0.0,0.0),(1.0, 0.25, 0.25))}
	rgbwcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbwdict,256)
	rgbkcmap = mplt.colors.LinearSegmentedColormap('my_colormap',rgbkdict,256)
	
	
	SpikesCmap.set_bad(bgColor)
	cbar1 =mplt.colorbar.ColorbarBase(ax=ax1,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar1.set_label('Excitability (decreasing)',fontsize='small')
	cbar1.set_ticks(np.array([0,0.5,1.0]))
	cbar1.set_ticklabels(np.array([0,0.5,1]))
	#cbar1.set_ticklabels(np.array([cMin,(cMin+cMax)/2,cMax]))
	cbar2 =mplt.colorbar.ColorbarBase(ax=ax3,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar2.set_label('Excitability (increasing)',fontsize='small')
	cbar2.set_ticks(np.array([0,0.5,1.0]))
	cbar2.set_ticklabels(np.array([0,0.5,1]))
	#ax.imshow((Q*1./(QN+(QN==0))[:,None]).transpose()\
	#ax.imshow(S.transpose()\
	#,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	#,interpolation='none',origin='lower',extent=(-100,100,0,100))
	#ax.set_ylim(0,100)
	#ax.set_xlim(-100,100)
	#ax.set_xticks(np.array([-100,0,100]))
	#ax.set_yticks(np.array([0,50,100]))np.reshape(W,((yMax-yMin)*Res,-1,3))
	ax.imshow(Y1, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax.set_ylim(yMin,yMax)
	ax.set_xlim(xMin,xMax)
	ax.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	#second plot (fraction considered)
	ax2.imshow(Y\
	,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.set_ylim(yMin,yMax)
	ax2.set_xlim(xMin,xMax)
	ax2.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	ax2.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	pylab.savefig(FileName)
	pylab.close()
	c=np.histogram2d(np.sign(b)*50+a,np.sign(b1)*50+a1,bins=(np.arange(-50,150),np.arange(-50,150)))[0]
	fig2=pylab.figure(1,figsize=(8,8))
	fig2.text(0.05,0.5,'typical excitability rank',rotation='vertical',ha='center',va='center')#,fontsize='xx-large'
	fig2.text(0.5,0.05,'typical excitability change',ha='center',va='center')#,fontsize='xx-large'
	fig2.text(0.97,0.3,'neg.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	fig2.text(0.97,0.7,'pos.deviation',rotation='vertical',fontsize='small',ha='center',va='center')
	fig2.text(0.3,0.97,'neg. change',ha='center',va='center',fontsize='small')
	fig2.text(0.7,0.97,'pos. change',ha='center',va='center',fontsize='small')
	ax4=fig2.add_axes([0.1, 0.1, 0.8, 0.8],axisbg=bgColor)
	ax5=fig2.add_axes([0.1, 0.9, 0.4, 0.03])
	ax6=fig2.add_axes([0.5, 0.9, 0.4, 0.03])
	ax7=fig2.add_axes([0.9, 0.1, 0.03, 0.4])
	ax8=fig2.add_axes([0.9, 0.5, 0.03, 0.4])
	ax4.set_xticks(np.array([]))
	ax4.set_yticks(np.array([]))
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	ax6.set_xticks(np.array([]))
	ax6.set_yticks(np.array([]))
	ax7.set_xticks(np.array([]))
	ax7.set_yticks(np.array([]))
	ax8.set_xticks(np.array([]))
	ax8.set_yticks(np.array([]))
	cbar5 =mplt.colorbar.ColorbarBase(ax=ax5,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar6 =mplt.colorbar.ColorbarBase(ax=ax6,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.), orientation='horizontal')
	cbar7 =mplt.colorbar.ColorbarBase(ax=ax7,cmap=rgbwcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar8 =mplt.colorbar.ColorbarBase(ax=ax8,cmap=rgbkcmap,norm=mplt.colors.Normalize(vmin=0,vmax=1.))
	cbar5.set_ticks(np.array([]))
	cbar6.set_ticks(np.array([]))
	cbar7.set_ticks(np.array([]))
	cbar8.set_ticks(np.array([]))
	ax4.imshow(np.log10(np.clip(c.transpose(),1,1e3)), aspect='equal', vmin=0., vmax=3., cmap=SpikesCmap\
	, interpolation='none', origin='lower', extent=(0,200,0,200))
	ax4.set_ylim(0,200)
	ax4.set_xlim(0,200)
	#ax4.set_xticks(np.array([xMin,(xMin+xMax)/2,xMax]))
	#ax4.set_yticks(np.array([yMin,(yMin+yMax)/2,yMax]))
	ax4.plot(np.array([100,100]),np.array([0,200]),'k-')
	ax4.plot(np.array([0,200]),np.array([100,100]),'k-')
	pylab.savefig(FileName2)
	pylab.close()
	#pylab.figure()
	#pylab.imshow(np.reshape(Hsum,((yMax-yMin)*Res,-1))\
	#,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	#pylab.colorbar()
	#pylab.show()
	return


def TransitionMatricesPlotSlices(HdfFile1, HdfFile2, FileName\
, FieldName='Excitability3/relE'\
, IndFieldName='Excitability3/Indices', LocFieldName='Locations', AmpFieldName='Amplitudes'\
, xMin=0, xMax=64, yMin=0, yMax=64, Res=8, cLabel=''\
, ampThreshold=2.5, cThreshold=-1e12\
, bgColor=pylab.cm.gray(0.5), SpikesCmap=pylab.cm.hot_r):
	#weighting? -- number of spikes? channels above 0.1 Hz?
	#representation as a (coloured:firingrate) scatterplot? or s.th continuous?
	#estimate a diffusion constant and a drift term per location?
	#compute histogram intersection per location? -- meaning??? stats???
	assert xMin<xMax
	#assert Res>=2
	A=np.array(['y-position','x-position'])
	fig=pylab.figure(1,figsize=(12,12))
	ax0=fig.add_axes([0.02, 0.02, 0.3, 0.3])
	ax1=fig.add_axes([0.35, 0.02, 0.3, 0.3])
	ax2=fig.add_axes([0.68, 0.02, 0.3, 0.3])
	ax3=fig.add_axes([0.02, 0.35, 0.3, 0.3])
	ax4=fig.add_axes([0.35, 0.35, 0.3, 0.3])
	ax5=fig.add_axes([0.68, 0.35, 0.3, 0.3])
	ax6=fig.add_axes([0.02, 0.68, 0.3, 0.3])
	ax7=fig.add_axes([0.35, 0.68, 0.3, 0.3])
	ax8=fig.add_axes([0.68, 0.68, 0.3, 0.3])
	ax0.set_xticks(np.array([]))
	ax0.set_yticks(np.array([]))
	ax1.set_xticks(np.array([]))
	ax1.set_yticks(np.array([]))
	ax2.set_xticks(np.array([]))
	ax2.set_yticks(np.array([]))
	ax3.set_xticks(np.array([]))
	ax3.set_yticks(np.array([]))
	ax4.set_xticks(np.array([]))
	ax4.set_yticks(np.array([]))
	ax5.set_xticks(np.array([]))
	ax5.set_yticks(np.array([]))
	ax6.set_xticks(np.array([]))
	ax6.set_yticks(np.array([]))
	ax7.set_xticks(np.array([]))
	ax7.set_yticks(np.array([]))
	ax8.set_xticks(np.array([]))
	ax8.set_yticks(np.array([]))
	g=h5py.File(HdfFile1,'r')
	Sampling=g['Sampling'].value
	Ind1=g[IndFieldName].value
	Amp1=g[AmpFieldName].value
	Amp1=Amp1[Ind1]
	cQ1=g[FieldName].value
	cQ1=cQ1[Ind1]
	pAmpInd1=(Amp1>=ampThreshold)*(cQ1>=cThreshold)
	cQ1=cQ1[pAmpInd1]
	Loc1=g[LocFieldName].value[Ind1,:][pAmpInd1,:]
	Amp1=Amp1[pAmpInd1]
	g.close()
	g=h5py.File(HdfFile2,'r')
	Ind2=g[IndFieldName].value
	Amp2=g[AmpFieldName].value
	Amp2=Amp2[Ind2]
	cQ2=g[FieldName].value
	cQ2=cQ2[Ind2]
	pAmpInd2=(Amp2>=ampThreshold)*(cQ2>=cThreshold)
	cQ2=cQ2[pAmpInd2]
	Loc2=g[LocFieldName].value[Ind2,:][pAmpInd2,:]
	Amp2=Amp2[pAmpInd2]
	g.close()
	#Excitability measure
	A1=np.argsort(np.argsort(100.*np.exp(-cQ1)))*100./len(cQ1)
	A2=np.argsort(np.argsort(100.*np.exp(-cQ2)))*100./len(cQ2)
	L1=len(A1)
	L2=len(A2)
	#histogram of excitability of different locations
	H1=np.histogramdd(np.concatenate((Loc1,A1[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H2=np.histogramdd(np.concatenate((Loc2,A2[:,None]),axis=1)\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res\
	,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res,np.arange(101)))[0]
	H1=H1[1:,1:,:]+H1[:-1,1:,:]+H1[:-1,:-1,:]+H1[1:,:-1,:]
	H2=H2[1:,1:,:]+H2[:-1,1:,:]+H2[:-1,:-1,:]+H2[1:,:-1,:]
	H1r=np.reshape(H1,(-1,100))
	H2r=np.reshape(H2,(-1,100))
	#histogram of spike count
	HN1=np.histogram2d(Loc1[:,0],Loc1[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN2=np.histogram2d(Loc2[:,0],Loc2[:,1]\
	,bins=((np.arange(yMin*Res,yMax*Res+2)-0.5)*1./Res,(np.arange(xMin*Res,xMax*Res+2)-0.5)*1./Res))[0]
	HN1=HN1[1:,1:]+HN1[:-1,1:]+HN1[:-1,:-1]+HN1[1:,:-1]
	HN2=HN2[1:,1:]+HN2[:-1,1:]+HN2[:-1,:-1]+HN2[1:,:-1]
	HN1=HN1.flatten()
	HN2=HN2.flatten()
	#threshold spike count
	#NInd=(HN1+HN2)>=100
	NInd=np.argsort(HN1+HN2)[-500:]
	Hsum=np.clip(np.sum((H2r-H1r)*np.arange(100)[None,:],axis=1)*1./np.clip((HN2+HN1),1,1e12),-100,100)
	##Hloc=np.histogram(Hsum,bins=np.arange(-100,101))[0]
	#Histogram of average and single location activity change
	print Hsum.shape
	Q=np.histogram2d((Hsum[:,None]*np.ones(100)[None,:]).flatten()\
	,(np.ones(Hsum.shape[0])[:,None]*np.arange(100)[None,:]).flatten()\
	,bins=(np.arange(-100,101),np.arange(101)),weights=(H2r-H1r).flatten())[0]
	QN=np.histogram(Hsum,bins=(np.arange(-100,101)),weights=np.sqrt(np.sum(H2r+H1r,axis=1)))[0]#
	#need to figure out where the largest difference is (in Ex.) take 10% window (sliding, periodic boundary)
	#use absolute difference (^1/3?)
	#cannot plot that...clustering? -- distances? - correlations... or firing rate cutoff take highest 500?
	#S=((H2r-H1r)*1./np.clip((np.sqrt(np.sum(H2r+H1r,axis=1))),1,1e12)[:,None])[NInd,:][np.argsort(Hsum[NInd]),:]
	S=np.concatenate((np.zeros((Hsum.shape[0],1)),np.cumsum(((H2r-H1r)\
	*1./np.clip((np.sqrt(np.sum(H2r+H1r,axis=1))),1,1e12)[:,None]),axis=1)),axis=1)
	Sinc=(S[:,20::10]-S[:,:-20][:,::10])/10.
	ax0.imshow(np.reshape(Sinc[:,0],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax1.imshow(np.reshape(Sinc[:,1],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax2.imshow(np.reshape(Sinc[:,2],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax3.imshow(np.reshape(Sinc[:,3],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax4.imshow(np.reshape(Sinc[:,4],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax5.imshow(np.reshape(Sinc[:,5],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax6.imshow(np.reshape(Sinc[:,6],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax7.imshow(np.reshape(Sinc[:,7],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	ax8.imshow(np.reshape(Sinc[:,8],(yMax*Res-yMin*Res,xMax*Res-xMin*Res))\
	,vmin=-1., vmax=1.,cmap=SpikesCmap, aspect='equal'\
	,interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	pylab.savefig(FileName)
	pylab.close()
	#pylab.figure()
	#pylab.imshow(np.reshape(Hsum,((yMax-yMin)*Res,-1))\
	#,aspect='equal',interpolation='none',origin='lower',extent=(xMin,xMax,yMin,yMax))
	#pylab.colorbar()
	#pylab.show()
	return
