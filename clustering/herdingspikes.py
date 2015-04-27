# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 11:17:38 2014

@author: Martino Sorbaro
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA
from scipy.stats import itemfreq
import h5py
import warnings
from sys import stdout


def ImportInterpolated(filename,shapesupto=None):
    g=h5py.File(filename,'r')
    A = spikeclass(np.array(g['Locations'].value,dtype=float).T)
    A.LoadTimes(g['Times'].value)
    A.SetSampling(g['Sampling'].value)
    if shapesupto == None:
        A.LoadShapes(np.array(g['Shapes'].value).T)
    elif shapesupto > 0:
        A.LoadShapes(np.array(g['Shapes'].value).T[:shapesupto])
    g.close()
    return A


class spikeclass(object):
    """A class containing code to work on 2d data with the Mean Shift clustering algorithms and various filters.

    Can be initialised in three ways:
        -- with a string pointing to an hdf5 file. This file must have been previously saved using this class.
        -- with a single [2,N] array containing raw data.
        -- with two arguments: a [2,Ndata] array containing data and a [Ndata] array containing, for every point, the ID of the cluster it belongs to.
        """

    def __init__(self,*args,**kwargs):
        if len(args) == 1:
            if isinstance(args[0],str):
                g=h5py.File(args[0],'r')
                self.__data=np.array(g['data'].value,dtype=float)
                self.__ClusterID = np.array(g['cluster_id'].value,dtype=int) if 'cluster_id' in g.keys() else np.array([])
                self.__c=np.array(g['centres'].value,dtype=float) if 'centres' in g.keys() else np.array([])
                self.__times = np.array(g['times']) if 'times' in g.keys() else np.array([])
                self.__shapes = np.array(g['shapes']) if 'shapes' in g.keys() else np.array([])
                self.__colours = np.array([])
                self.__sampling = g['Sampling'].value
                g.close()
            else:
                givendata = args[0]
                ndata = np.shape(givendata)[1]
                if np.shape(givendata) != (2,ndata): raise ValueError('Data must be a (2,N) array')
                self.__data = givendata
                self.__ClusterID = np.array([])
                self.__c = []
                self.__times = np.array([])
                self.__shapes = np.array([])
                self.__colours = np.array([])
                self.__sampling = []
        elif len(args) == 2:
            ndata = args[0].shape[1]
            if np.shape(args[0]) != (2,ndata): raise ValueError('Data must be a (2,N) array')
            self.__data = args[0]
            self.__c = np.zeros([2,np.max(args[1])+1])
            self.__ClusterID = np.array(args[1])
            self.__times = np.array([])
            self.__shapes = np.array([])
            self.__colours = np.array([])
            self.__sampling = []
        else:
            raise ValueError('Can be initialised with 1 argument (the data set or a file) or 2 arguments (data, ClusterID)')
        self.Backup()


    def Backup(self):
        """Creates a checkpoint, to be used for a subsequent call to UndoLast()"""
        self.__backup = {0:self.__data,1:self.__ClusterID,2:self.__c,3:self.__shapes,4:self.__times}


    def UndoLast(self):
        """The object restores the data as it was before the last call of a filter, or Backup()."""
        self.__data,self.__ClusterID,self.__c,self.__shapes,self.__times = self.__backup[0],self.__backup[1],self.__backup[2],self.__backup[3],self.__backup[4]

    def Colours(self):
        if np.shape(self.__colours)[0] != self.NClusters():
            colours=plt.cm.spectral(np.random.permutation(np.linspace(0,1,num=self.NClusters())))
            colours=np.append(np.array([0,0,0,0.5]),colours[:-1])
            self.__colours=np.reshape(colours,(self.NClusters(),4))
        return self.__colours

    def LogHistPlot(self,save=None):
        """Plots a density histogram. This does not work with the current data, but with data loaded by the most recent class constructor call."""
        #fig,ax = plt.subplots()
        #ax.imshow(np.log10(self.__hist),origin='lower',interpolation='none')
        #plt.show()
        plt.figure(figsize=(8,8))
        ax = plt.subplot(111)
        ax.set_axis_bgcolor('black')
        n,xb,yb=np.histogram2d(self.__data[1], self.__data[0], bins=(np.arange(0,65,0.2),np.arange(0,65,0.2)))
        rateMasked = np.ma.array (n, mask=(n==0))
        cmap = plt.cm.RdBu_r
        cmap.set_bad('k')
        plt.pcolor(yb,xb,np.log(rateMasked),cmap=cmap)
        plt.axis('equal')
        plt.xlim((0,65))
        plt.ylim((0,65))
        if save != None:
            plt.savefig(save)

    def DataPlot(self):
        """Plots the current data. If clustering was performed, the cluster centres and ID (colour) are plotted, otherwise a black and white scatterplot is plotted"""
        fig,ax = plt.subplots(figsize=(7,7))
        ax.set_axis_bgcolor('black')
        if np.size(self.__ClusterID):
            ax.scatter(self.__data[0],self.__data[1],c=self.Colours()[self.__ClusterID],marker='o',s=1,edgecolors='none',alpha=0.5)
        else:
            ax.scatter(self.__data[0],self.__data[1],marker=',',c='w',s=1,edgecolors='none',alpha=0.5)
        ax.set_aspect('equal')
        ax.set_xlim([min(self.__data[0]),max(self.__data[0])])
        ax.set_ylim([min(self.__data[1]),max(self.__data[1])])

    def PartPlot(self,window,save=None):
        PartPlot(self,window[0],window[1],window[2],window[3],save)

    def PartPlot(self,x1,x2,y1,y2,save=None):
        ratio = (x2-x1)/(y2-y1)
        plt.figure(figsize=(12*ratio,12))
        ax = plt.subplot(121)
        ax.grid(which='major', axis='x', linewidth=1, linestyle='-', color='0.75')
        ax.grid(which='major', axis='y', linewidth=1, linestyle='-', color='0.75')
        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)
        plt.scatter(self.__data[0],self.__data[1],marker='o',s=3,edgecolors='none',c=self.Colours()[self.__ClusterID])
        ax.set_aspect('equal')
        plt.xticks(np.arange(np.round(x1),np.ceil(x2)))
        plt.yticks(np.arange(np.round(y1),np.ceil(y2)))
        if save != None:
            plt.savefig(save)

    def ShapesPlot(self,clusters=None,save=None):
        if clusters == None:
            clusters = range(self.NClusters())
        if len(clusters)>12:
            warnings.warn("Only the first 12 of the given clusters are shown")
            clusters = clusters[:12]
        plt.figure(1, figsize=(20, 2.5))
        ax=plt.subplot(111, frameon=False)
        ax.grid(which='major', axis='x', linewidth=1, linestyle='-', color='0.75')
        ax.grid(which='major', axis='y', linewidth=1, linestyle='-', color='0.75')

        clShapes = lambda c: self.__shapes[:,self.__ClusterID == c]
        sl = 1.0*np.shape(self.__shapes)[0]

        for ic,c in enumerate(clusters):
            myShapes = clShapes(c)
            plInds=range(np.min([30,myShapes.shape[1]]))
            [plt.plot(ic+np.arange(sl)/sl-.5,myShapes[:,i],color=self.Colours()[c],alpha=0.2) for i in plInds]
            plt.plot(ic+np.arange(sl)/sl-.5,np.mean(myShapes,axis=1),'-',color=self.__colours[c],lw=2.5)

        plt.xlim((-.5,11.5))
        plt.yticks([])
        plt.xticks(np.arange(np.min((len(clusters),12))))
        plt.xlabel('Cluster ID')
        plt.grid(0)
        if save != None:
            plt.savefig(save)

    def NData(self):
        """Returns the current number of datapoints."""
        return np.shape(self.__data)[1]

    def NClusters(self):
        """Returns the current number of clusters, or 0 if no clustering was performed."""
        if np.size(self.__ClusterID):
            return np.shape(self.__c)[1]
        else:
            return 0

    def Locations(self):
        """Returns the data set."""
        return self.__data

    def Shapes(self):
        """Returns the shapes set."""
        return self.__shapes

    def Times(self):
        """Returns the times set."""
        return self.__times

    def ClusterID(self):
        """Returns an array containing the id of the cluster every data point belongs to."""
        return self.__ClusterID

    def ClusterLoc(self):
        """Returns an array containing the locations of the cluster centres."""
        return np.array(self.__c)

    def ClusterSizes(self):
        """Returns an array containing the number of points in each cluster."""
        return itemfreq(self.__ClusterID)[:,1]

    def Sampling(self):
        """Returns the sampling rate."""
        return self.__sampling

    def LoadShapes(self, shapes):
        assert np.size(np.shape(shapes)) == 2
        assert np.shape(shapes)[1] == self.NData()
        self.__shapes = np.array(shapes)

    def LoadTimes(self, times):
        assert np.size(np.shape(times)) == 1
        assert np.shape(times)[0] == self.NData()
        self.__times = np.array(times)

    def SetSampling(self, s):
        self.__sampling = s

    def Save(self,string):
        """Saves data, cluster centres and ClusterIDs to a hdf5 file."""
        g=h5py.File(string+'.hdf5','w-')
        g.create_dataset("data",data=self.__data)
        if self.__c != np.array([]):
            g.create_dataset("centres",data=self.__c)
        if self.__ClusterID != np.array([]):
            g.create_dataset("cluster_id",data=self.__ClusterID)
        if self.__times != np.array([]):
            g.create_dataset("times",data=self.__times)
        if self.__shapes != np.array([]):
            g.create_dataset("shapes",data=self.__shapes)
        g.close()

    def MeanShift(self,h):
        """Performs the scikit-learn Mean Shift clustering. kwargs are passed to the MeanShift class."""
        MS = MeanShift(bin_seeding=True,bandwidth=h, cluster_all=True, min_bin_freq=0)
        print "Starting sklearn Mean Shift... ",
        stdout.flush()
        MS.fit_predict(self.__data.T)
        self.__ClusterID = MS.labels_
        self.__c = MS.cluster_centers_.T
        print "done."
        stdout.flush()

    def ShapePCA(self,ncomp=None,white=False):
        print "Starting sklearn PCA...",
        stdout.flush()
        p = PCA(n_components=ncomp,whiten=white)
        fit = p.fit_transform(self.Shapes().T).T
        print "done."
        stdout.flush()
        return fit

    def CombinedMeanShift(self,h,alpha,PrincComp=[]):
        MS = MeanShift(bin_seeding=True, bandwidth=h, cluster_all=True, min_bin_freq=0)
        if not PrincComp.any():
          PrincComp = self.ShapePCA(1)
        print "Starting sklearn Mean Shift... ",
        stdout.flush()
        print self.__data.shape, PrincComp.shape
        fourvector = np.vstack((self.__data,alpha*PrincComp))
        MS.fit_predict(fourvector.T)
        self.__ClusterID = MS.labels_
        self.__c = MS.cluster_centers_.T[0:1]
        print "done."
        stdout.flush()

    def RemoveData(self,newn):
        """Randomly chooses datapoints and deletes all the others

        Arguments:

        newn -- the number of datapoints to be kept"""
        self.Backup()
        initialn = self.NData()
        if newn<self.NData():
            ind = np.random.choice(range(self.NData()),size=newn,replace=False)
            self.KeepOnly(ind)
            print('RemoveData removed '+str(initialn-self.NData())+' datapoints.')
        else:
            print('RemoveData: No points were discarded.')

    # A Dangerous Method. Do not use more than once.
    # Also, I think it's better not to use it in combination with ReduceData
    def FilterLowDensity(self,threshold,nbins=[400,400]):
        """Bins points in 100 bins per axis and deletes points in bins with number of points <= threshold.

        Returns an array containing the indices corresponding to KEPT data.
        """
        self.Backup()
        hist,bx,by = np.histogram2d(self.__data[0],self.__data[1],nbins)
        # the *1.001 is needed to include the rightmost and topmost points in the bins ... bad coding indeed.
        binspanx = (np.max(self.__data[0])-np.min(self.__data[0]))/nbins[0]*1.001
        binspany = (np.max(self.__data[1])-np.min(self.__data[1]))/nbins[1]*1.001
        nbx = ((self.__data[0]-np.min(self.__data[0]))//binspanx).astype(int)
        nby = ((self.__data[1]-np.min(self.__data[1]))//binspany).astype(int)
        initialn = self.NData()
        ind = np.where(hist[nbx,nby]>threshold)[0]
        self.KeepOnly(ind)
        print('FilterLowDensity removed '+str(initialn-self.NData())+' datapoints.')
        return ind


    def FilterSmallClusters(self,threshold):
        """Removes all datapoints belonging to clusters with 'threshold' or less datapoints."""
        self.Backup()
        numclus = self.NClusters()
        initialdata = self.NData()
        sizes = self.ClusterSizes()
        #create a conversion table to get rid of gaps in cluster IDs
        c_ind_kept = np.where(sizes >= threshold)[0]
        newID = -np.ones(numclus,dtype=np.int)
        newID[c_ind_kept] = np.array(range(len(c_ind_kept)))
        #update temporarily the ClusterID vector
        self.__ClusterID = newID[self.__ClusterID]
        #delete data whose cluster was deleted, and clusters
        d_ind_kept = np.where(self.__ClusterID != -1)[0]
        self.KeepOnly(d_ind_kept)
        self.__ClusterID = self.__ClusterID[d_ind_kept]
        self.__c = self.__c[:,c_ind_kept]
        print('FilterSmallClusters removed '+str(numclus-self.NClusters())+' clusters and '+str(initialdata-self.NData())+' datapoints.')
        return d_ind_kept

    def KeepOnly(self,ind_kept):
        #does not act on clusters!
        self.__data = self.__data[:,ind_kept]
        if np.size(self.__shapes):
            self.__shapes = self.__shapes[:,ind_kept]
        if np.size(self.__times):
            self.__times = self.__times[ind_kept]

    def CropClusters(self,xmin,xmax,ymin,ymax,keepinside=True):
        """Removes all datapoints belonging to clusters centered outside the specified area."""
        self.Backup()
        numclus = self.NClusters()
        initialdata = self.NData()
        cx,cy = self.__c
        #create a conversion table to get rid of gaps in cluster IDs
        if keepinside:
            condition = [ x&y&z&w for (x,y,z,w) in zip(cx<=xmax,cx>=xmin,cy<=ymax,cy>=ymin)]
        else:
            condition = [ -(x&y&z&w) for (x,y,z,w) in zip(cx<=xmax,cx>=xmin,cy<=ymax,cy>=ymin)]
        c_ind_kept = np.where(condition)[0]
        newID = -np.ones(numclus,dtype=np.int)
        newID[c_ind_kept] = np.array(range(len(c_ind_kept)))
        #update temporarily the ClusterID vector
        self.__ClusterID = newID[self.__ClusterID]
        #delete data whose cluster was deleted, and clusters
        d_ind_kept = np.where(self.__ClusterID != -1)[0]
        self.KeepOnly(d_ind_kept)
        self.__ClusterID = self.__ClusterID[d_ind_kept]
        self.__c = self.__c[:,c_ind_kept]
        print('CropClusters removed '+str(numclus-self.NClusters())+' clusters and '+str(initialdata-self.NData())+' datapoints.')
        return d_ind_kept

    def Crop(self,xmin,xmax,ymin,ymax,keepinside=True):
        # TO BE TESTED
        self.Backup()
        dx,dy = self.__data
        numclus = self.NClusters()
        initialdata = self.NData()
        if keepinside:
            condition = [ x&y&z&w for (x,y,z,w) in zip(dx<=xmax,dx>=xmin,dy<=ymax,dy>=ymin)]
        else:
            condition = [ -x&y&z&w for (x,y,z,w) in zip(dx<=xmax,dx>=xmin,dy<=ymax,dy>=ymin)]
        d_ind_kept = np.where(condition)[0]
        if numclus:
            cid_kept_all = self.__ClusterID[d_ind_kept]
            c_ind_kept = np.unique(cid_kept_all)
            newID = -np.ones(numclus,dtype=np.int)
            newID[c_ind_kept] = np.array(range(len(c_ind_kept)))
            #update temporarily the ClusterID vector
            self.__ClusterID = newID[self.__ClusterID]
            #delete data whose cluster was deleted, and clusters
            self.__ClusterID = self.__ClusterID[d_ind_kept]
            self.__c = self.__c[:,c_ind_kept]
        self.KeepOnly(d_ind_kept)
        print('Crop removed '+str(numclus-self.NClusters())+' clusters and '+str(initialdata-self.NData())+' datapoints.')
        return d_ind_kept

    def Classify(self,nbins = [40,40],threshold = 5):
        hist,bx,by = np.histogram2d(self.__data[0],self.__data[1],nbins)
        binspanx = (np.max(self.__data[0])-np.min(self.__data[0]))/nbins[0]*1.001
        binspany = (np.max(self.__data[1])-np.min(self.__data[1]))/nbins[1]*1.001
        nbx = ((self.__data[0]-np.min(self.__data[0]))//binspanx).astype(int)
        nby = ((self.__data[1]-np.min(self.__data[1]))//binspany).astype(int)
        ind = np.where(hist[nbx,nby]<=threshold)[0]
        print "Based on "+str(len(ind))+" examples of bad shapes."
        normalise = lambda X: X/np.max(np.abs(X),axis=0)
        badshape = np.mean(normalise(self.Shapes()),axis=1)
        score = np.dot(badshape,self.Shapes())
        scorePCA = self.ShapePCA(ncomp=1)[0]
        return score,scorePCA








