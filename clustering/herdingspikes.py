# -*- coding: utf-8 -*-
"""
Created on Tue Sep 23 11:17:38 2014

@author: Martino Sorbaro
@author: Matthias Hennig
"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from sklearn import __version__ as skvers
from sklearn.cluster import MeanShift
from sklearn.decomposition import PCA
from sklearn import svm, mixture
from sklearn.metrics.pairwise import euclidean_distances
from scipy.stats import itemfreq
import h5py
import warnings
from sys import stdout

if float(skvers) < 0.17:
    raise Warning('Sklearn version >= 0.17 may be needed')


def ImportInterpolated(filename, shapesupto=None):
    """Helper function to read spike data from an hdf5 file."""
    g = h5py.File(filename, 'r')
    A = spikeclass(np.array(g['Locations'].value, dtype=float).T)
    A.LoadTimes(np.floor(g['Times'].value).astype(int))
    A.SetSampling(g['Sampling'].value)
    if shapesupto is None:
        A.LoadShapes(np.array(g['Shapes'].value).T)
    elif shapesupto > 0:
        A.LoadShapes(np.array(g['Shapes'].value).T[:shapesupto])
    g.close()
    A._spikeclass__expinds = np.array([0])
    return A


def ImportInterpolatedList(filenames, shapesupto=22):
    """ Helper function to read in spike data from a list of hdf5 files.
    Returns a class object which keeps track of the
    indices where each file begins.
    """
    loc = np.array([[], []], dtype=float)
    t = np.array([], dtype=int)
    sh = np.array([], dtype=int)
    inds = np.zeros(len(filenames)+1, dtype=int)
    s = np.zeros(len(filenames))
    for i, f in enumerate(filenames):
        g = h5py.File(f, 'r')
        print(f)
        loc = np.append(loc, g['Locations'].value.T, axis=1)
        inds[i] = len(t)  # store index of first spike
        t = np.append(t, np.floor(g['Times'].value).astype(int))
        s[i] = g['Sampling'].value
        # if shapesupto == None:
        #  sh = np.append(sh, np.array(g['Shapes'].value))
        # elif shapesupto > 0:
        sh = np.append(sh, np.array(g['Shapes'].value)[:, :shapesupto])
        g.close()
    inds[len(filenames)] = len(t)
    sh = np.reshape(sh, (len(t), shapesupto))
    if len(np.unique(s)) > 1:
        raise Warning('Data sets have different sampling rates\n' + str(s))
    A = spikeclass(loc)
    A.LoadTimes(t)
    A.SetSampling(s[0])
    A.LoadShapes(sh.T)
    A._spikeclass__expinds = inds
    return A


class spikeclass(object):
    """A class containing code to work on 2d data with the Mean Shift
    clustering algorithms and various filters.

    Can be initialised in three ways:
        -- with a string pointing to an hdf5 file. This file must have been
         previously saved using this class.
        -- with a single [2,N] array containing raw data.
        -- with two arguments: a [2,Ndata] array containing data and a [Ndata]
         array containing, for every point,
         the ID of the cluster it belongs to.
        """

    def __init__(self, *args, **kwargs):
        if len(args) == 1:
            if isinstance(args[0], str):
                g = h5py.File(args[0], 'r')
                self.__data = np.array(g['data'].value, dtype=float)
                self.__ClusterID = np.array(g['cluster_id'].value, dtype=int) \
                    if 'cluster_id' in g.keys() else np.array([])
                self.__c = np.array(g['centres'].value, dtype=float) \
                    if 'centres' in g.keys() else np.array([])
                self.__times = np.array(g['times']) \
                    if 'times' in g.keys() else np.array([])
                self.__shapes = np.array(g['shapes']) \
                    if 'shapes' in g.keys() else np.array([])
                self.__colours = np.array([])
                self.__sampling = g['Sampling'].value \
                    if 'Sampling' in g.keys() else np.array([])
                self.__expinds = g['expinds'].value \
                    if 'expinds' in g.keys() else np.array([0])
                g.close()
            else:
                givendata = args[0]
                ndata = np.shape(givendata)[1]
                if np.shape(givendata) != (2, ndata):
                    raise ValueError('Data must be a (2,N) array')
                self.__data = givendata
                self.__ClusterID = np.array([])
                self.__c = []
                self.__times = np.array([])
                self.__shapes = np.array([])
                self.__colours = np.array([])
                self.__sampling = []
                self.__expinds = np.array([0])
        elif len(args) == 2:
            ndata = args[0].shape[1]
            if np.shape(args[0]) != (2, ndata):
                raise ValueError('Data must be a (2,N) array')
            self.__data = args[0]
            self.__c = np.zeros([2, np.max(args[1])+1])
            self.__ClusterID = np.array(args[1])
            self.__times = np.array([])
            self.__shapes = np.array([])
            self.__colours = np.array([])
            self.__sampling = []
            self.__expinds = np.array([0])
        else:
            raise ValueError(
                'Can be initialised with 1 argument (the data' +
                ' set or a file) or 2 arguments (data, ClusterID)')
        self.Backup()

    def Colours(self):
        if np.shape(self.__colours)[0] != self.NClusters():
            colours = plt.cm.spectral(np.random.permutation(np.linspace(
                0, 1, num=self.NClusters())))
            colours = np.append(np.array([0, 0, 0, 0.5]), colours[:-1])
            self.__colours = np.reshape(colours, (self.NClusters(), 4))
        return self.__colours

# PLOTTING METHODS

    def LogHistPlot(self, save=None):
        """Plots a density histogram."""
        # fig,ax = plt.subplots()
        # ax.imshow(np.log10(self.__hist),origin='lower',interpolation='none')
        # plt.show()
        plt.figure(figsize=(8, 8))
        ax = plt.subplot(111)
        ax.set_axis_bgcolor('black')
        n, xb, yb = np.histogram2d(
            self.__data[1], self.__data[0],
            bins=(np.arange(0, 65, 0.2), np.arange(0, 65, 0.2)))
        rateMasked = np.ma.array(n, mask=(n == 0))
        cmap = plt.cm.RdBu_r
        cmap.set_bad('k')
        plt.pcolor(yb, xb, np.log(rateMasked), cmap=cmap)
        plt.axis('equal')
        plt.xlim((0, 65))
        plt.ylim((0, 65))
        if save is not None:
            plt.savefig(save)

    def DataPlot(self, save=None):
        """Plots the current data. If clustering was performed,
         the cluster centres and ID (colour) are plotted,
         otherwise a black and white scatterplot is plotted"""
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.set_axis_bgcolor('black')
        if np.size(self.__ClusterID):
            ax.scatter(self.__data[0], self.__data[1],
                       c=self.Colours()[self.__ClusterID], marker='o',
                       s=1, edgecolors='none', alpha=0.5)
        else:
            ax.scatter(self.__data[0], self.__data[1], marker=',',
                       c='w', s=1, edgecolors='none', alpha=0.5)
        ax.set_aspect('equal')
        ax.set_xlim([min(self.__data[0]), max(self.__data[0])])
        ax.set_ylim([min(self.__data[1]), max(self.__data[1])])
        if save is not None:
            plt.savefig(save)

    def PartPlot(self, rectangle, save=None):
        """Plots a portion of the space. Doesn't work prior to clustering."""
        (x1, x2, y1, y2) = rectangle
        ratio = (x2-x1)/(y2-y1)
        plt.figure(figsize=(12*ratio, 12))
        ax = plt.subplot(121)
        ax.grid(which='major', axis='x', linewidth=1,
                linestyle='-', color='0.75')
        ax.grid(which='major', axis='y', linewidth=1,
                linestyle='-', color='0.75')
        ax.set_xlim(x1, x2)
        ax.set_ylim(y1, y2)
        plt.scatter(self.__data[0], self.__data[1], marker='o', s=3,
                    edgecolors='none',
                    c=self.Colours()[self.__ClusterID])
        ax.set_aspect('equal')
        plt.xticks(np.arange(np.round(x1), np.ceil(x2)))
        plt.yticks(np.arange(np.round(y1), np.ceil(y2)))
        if save is not None:
            plt.savefig(save)

    def SpikesInCluster(self, c):
        return self.__ClusterID == c

    def ShapesPlot(self, clusters=None, save=None):
        """Plots the shapes of up to 12 clusters."""
        if clusters is None:
            clusters = range(self.NClusters())
        if len(clusters) > 12:
            warnings.warn("Only the first 12 of the given clusters are shown")
            clusters = clusters[:12]
        plt.figure(1, figsize=(20, 2.5))
        ax = plt.subplot(111, frameon=False)
        ax.grid(which='major', axis='x', linewidth=1,
                linestyle='-', color='0.75')
        ax.grid(which='major', axis='y', linewidth=1,
                linestyle='-', color='0.75')

        sl = 1.0*np.shape(self.__shapes)[0]

        for ic, c in enumerate(clusters):
            myShapes = self.__shapes[:, self.SpikesInCluster(c)]
            plInds = range(np.min([30, myShapes.shape[1]]))
            [plt.plot(ic+np.arange(sl)/sl-.5, myShapes[:, i],
             color=self.Colours()[c], alpha=0.2) for i in plInds]
            plt.plot(ic+np.arange(sl)/sl-.5, np.mean(myShapes, axis=1),
                     '-', color=self.__colours[c], lw=2.5)

        plt.xlim((-.5, 11.5))
        plt.yticks([])
        plt.xticks(np.arange(np.min((len(clusters), 12))))
        ax.set_xticklabels(clusters)
        plt.xlabel('Cluster ID')
        plt.grid(0)
        if save is not None:
            plt.savefig(save)

# GET AND SET METHODS

    def NData(self):
        """Returns the current number of datapoints."""
        return np.shape(self.__data)[1]

    def NClusters(self):
        """Returns the current number of clusters,
        or 0 if no clustering was performed."""
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
        """Returns an array containing the id of the cluster
        every data point belongs to."""
        return self.__ClusterID

    def ClusterLoc(self):
        """Returns an array containing the locations of the cluster centres."""
        return np.array(self.__c)

    def ClusterSizes(self):
        """Returns an array containing the number of points in each cluster."""
        return itemfreq(self.__ClusterID)[:, 1]

    def Sampling(self):
        """Returns the sampling rate."""
        return self.__sampling

    def LoadShapes(self, shapes):
        """Loads a KxN array, where K is the length of a single wave
        and N is the number of spikes, in the shapes vector."""
        assert np.size(np.shape(shapes)) == 2
        assert np.shape(shapes)[1] == self.NData()
        self.__shapes = np.array(shapes)

    def LoadTimes(self, times):
        """Loads a vector of spike times."""
        assert np.size(np.shape(times)) == 1
        assert np.shape(times)[0] == self.NData()
        self.__times = np.array(times, dtype=int)

    def SetSampling(self, s):
        """Sets the value of the sampling rate for internal usage."""
        self.__sampling = s

    def ExperimentIndices(self, i):
        """Returns a pair of indices denoting the start and end
        of an experiment. Can currently only be used if data from multiple
        experiments is read with the helper function ImportInterpolatedList."""
        raise NotImplementedError()
        # if i+1 < len(self.__expinds):
        #     final = self.__expinds[i+1]
        # elif i+1 == len(self.__expinds):
        #     final = self.NData()
        # else:
        #     raise ValueError('There are only ' + len(self.__expinds) +
        #                      ' datasets.')
        # return np.arange(self.__expinds[i], self.__expinds[i+1])

    def ClusterIndices(self, n, exper=None):
        raise NotImplementedError()
        # TO BE TESTED
        # idx = np.where(self.__ClusterID == n)[0]
        # if exper is not None:
        #     if exper + 1 < len(self.__indices):
        #         endind = self.__expinds[i+1]
        #         startind = self.__expinds[i]
        #     elif exper + 1 == len(self.__indices):
        #         endind = self.NData()
        #         startind = self.__expinds[i]
        #     else:
        #         raise ValueError('There are only ' + len(self.__indices) +
        #                          ' datasets.')
        #     idx = idx[idx >= startind]
        #     idx = idx[idx < endind]
        # return idx

    def ExperimentHeads(self):
        return self.__expinds

# SAVE

    def Save(self, string, compression=None):
        """Saves data, cluster centres and ClusterIDs to a hdf5 file.
        Offers compression of the shapes, 'lzf'
        appears a good trade-off between speed and performance.'"""
        g = h5py.File(string, 'w')
        g.create_dataset("data", data=self.__data)
        g.create_dataset("expinds", data=self.__expinds)
        if self.__c != np.array([]):
            g.create_dataset("centres", data=self.__c)
        if self.__ClusterID != np.array([]):
            g.create_dataset("cluster_id", data=self.__ClusterID)
        if self.__times != np.array([]):
            g.create_dataset("times", data=self.__times)
        if self.__shapes != np.array([]):
            g.create_dataset("shapes",
                             data=self.__shapes,
                             compression=compression)
        if self.__sampling:
            g.create_dataset("Sampling", data=self.__sampling)
        g.close()

# CLUSTERING AND ANALYSIS

    # def MeanShift(self, h, njobs=cpu_count()):
    #     """Performs the scikit-learn Mean Shift clustering.
    #     kwargs are passed to the MeanShift class."""
    #     MS = MeanShift(bin_seeding=True, bandwidth=h, cluster_all=True,
    #                    min_bin_freq=1, n_jobs=njobs)
    #     print("Starting sklearn Mean Shift... ")
    #     stdout.flush()
    #     MS.fit_predict(self.__data.T)
    #     self.__ClusterID = MS.labels_
    #     self.__c = MS.cluster_centers_.T
    #     print(self.__c)
    #     print("done.")
    #     stdout.flush()

    def AlignShapes(self):
        # Todo: optimise!
        """Re-aligns the peaks of the spike shapes. This can reduce spurious
        clustering at low sampling rates. Note the original shapes are
        overwritten and the resulting array is zero-padded at start and end.
        """
        peaks = np.argmin(self.Shapes(), axis=0)
        ap = int(np.median(peaks))
        peaks = -np.argmin(self.Shapes()[ap-2:ap+2], axis=0) + 1
        alShapes = np.insert(
            self.Shapes(),
            [0, 0, self.Shapes().shape[0], self.Shapes().shape[0]],
            0, axis=0)
        for d in np.arange(-2, 2):
            idxd = peaks == d
            alShapes[:, idxd] = np.roll(alShapes[:, idxd], d, axis=0)
        self.LoadShapes(alShapes)

    def ShapePCA(self, ncomp=None, white=False):
        """Compute PCA projections of spike shapes.
        If there are more than 1Mio data points, randomly sample 1Mio
        shapes and compute PCA from this subset only. Projections are
        then returned for all shapes.

        Arguments:

        ncomp -- the number of components to return
        white -- Perform whitening of data if set to True"""

        print("Starting sklearn PCA...")
        stdout.flush()
        p = PCA(n_components=ncomp, whiten=white)
        if self.NData() > 1000000:
            print(str(self.NData()) +
                  " spikes, using 1Mio shapes randomly sampled...")
            inds = np.random.choice(self.NData(), 1000000, replace=False)
            p.fit(self.Shapes()[:, inds].T)
            # compute projections
            fit = p.transform(self.Shapes().T).T
        else:
            print("using all "+str(self.NData())+" shapes...")
            fit = p.fit_transform(self.Shapes().T).T
        print("done.")
        stdout.flush()
        return fit

    def CombinedMeanShift(self, h, alpha,
                          PrincComp=None,
                          njobs=-2,
                          mbf=1):
        """Performs the scikit-learn Mean Shift clustering.

        Arguments:

        h -- the bandwidth
        alpha -- the weight of the principal components as compared
        to the spatial data.
        PrincComp -- used to pass already-computed principal components
        njobs -- the number of processes to be used (default: n. of CPU - 1)
        mbf -- the minimum number of items in a seed"""
        MS = MeanShift(bin_seeding=True, bandwidth=h, cluster_all=True,
                       min_bin_freq=mbf, n_jobs=njobs)
        if PrincComp is None:
            PrincComp = self.ShapePCA(2)
        print("Starting sklearn Mean Shift... ")
        stdout.flush()
        fourvector = np.vstack((self.__data, alpha*PrincComp))
        MS.fit_predict(fourvector.T)
        self.__ClusterID = MS.labels_
        self.__c = MS.cluster_centers_.T  # [0:1]
        print("done.")
        stdout.flush()

# FILTERS

    def RemoveData(self, newn):
        """Randomly chooses datapoints and deletes all the others

        Arguments:

        newn -- the number of datapoints to be kept"""
        self.Backup()
        initialn = self.NData()
        if newn < self.NData():
            ind = np.random.choice(range(self.NData()),
                                   size=newn, replace=False)
            self.KeepOnly(ind)
            print('RemoveData removed ' +
                  str(initialn-self.NData()) +
                  ' datapoints.')
        else:
            print('RemoveData: No points were discarded.')

    # A Dangerous Method. Do not use more than once.
    # Also, I think it's better not to use it in combination with ReduceData
    def FilterLowDensity(self, threshold, nbins=[400, 400]):
        """Bins points in 100 bins per axis and deletes points
        in bins with number of points <= threshold.

        Returns an array containing the indices corresponding to KEPT data.
        """
        self.Backup()
        hist, bx, by = np.histogram2d(self.__data[0], self.__data[1], nbins)
        # the *1.001 is needed to include the rightmost and topmost points
        # in the bins ... bad coding indeed.
        binspanx = (np.max(self.__data[0]) -
                    np.min(self.__data[0]))/nbins[0]*1.001
        binspany = (np.max(self.__data[1]) -
                    np.min(self.__data[1]))/nbins[1]*1.001
        nbx = ((self.__data[0]-np.min(self.__data[0]))//binspanx).astype(int)
        nby = ((self.__data[1]-np.min(self.__data[1]))//binspany).astype(int)
        initialn = self.NData()
        ind = np.where(hist[nbx, nby] > threshold)[0]
        self.KeepOnly(ind)
        print('FilterLowDensity removed ' +
              str(initialn-self.NData()) + ' datapoints.')
        return ind

    def FilterSmallClusters(self, threshold):
        """Removes all datapoints belonging to clusters with 'threshold'
        or less datapoints."""
        self.Backup()
        numclus = self.NClusters()
        initialdata = self.NData()
        sizes = self.ClusterSizes()
        # create a conversion table to get rid of gaps in cluster IDs
        c_ind_kept = np.where(sizes >= threshold)[0]
        newID = -np.ones(numclus, dtype=np.int)
        newID[c_ind_kept] = np.array(range(len(c_ind_kept)))
        # update temporarily the ClusterID vector
        self.__ClusterID = newID[self.__ClusterID]
        # delete data whose cluster was deleted, and clusters
        d_ind_kept = np.where(self.__ClusterID != -1)[0]
        self.KeepOnly(d_ind_kept)
        self.__ClusterID = self.__ClusterID[d_ind_kept]
        self.__c = self.__c[:, c_ind_kept]
        print('FilterSmallClusters removed '+str(numclus-self.NClusters()) +
              ' clusters and ' + str(initialdata-self.NData()) +
              ' datapoints.')
        return d_ind_kept

    def CropClusters(self, rectangle, outside=False, remove=True):
        """Keeps only datapoints belonging to clusters whose centres are
        inside the relevant window, or outside, if outside=True is passed.
        If remove=False, returns the IDs of the clusters in the area,
        without removing the rest."""
        (xmin, xmax, ymin, ymax) = rectangle
        numclus = self.NClusters()
        initialdata = self.NData()
        cx, cy = self.__c[:2]
        # create a conversion table to get rid of gaps in cluster IDs
        if not outside:
            condition = [x & y & z & w for (x, y, z, w) in
                         zip(cx <= xmax, cx >= xmin, cy <= ymax, cy >= ymin)]
        else:
            condition = [-(x & y & z & w) for (x, y, z, w) in
                         zip(cx <= xmax, cx >= xmin, cy <= ymax, cy >= ymin)]
        c_ind_kept = np.where(condition)[0]
        if remove:
            newID = -np.ones(numclus, dtype=np.int)
            newID[c_ind_kept] = np.array(range(len(c_ind_kept)))
            self.Backup()
            # update temporarily the ClusterID vector
            self.__ClusterID = newID[self.__ClusterID]
            # delete data whose cluster was deleted, and clusters
            d_ind_kept = np.where(self.__ClusterID != -1)[0]
            self.KeepOnly(d_ind_kept)
            self.__ClusterID = self.__ClusterID[d_ind_kept]
            self.__c = self.__c[:, c_ind_kept]
            print('CropClusters removed ' + str(numclus-self.NClusters()) +
                  ' clusters and ' + str(initialdata-self.NData()) +
                  ' datapoints.')
        return c_ind_kept

    def Crop(self, rectangle, outside=False, remove=True):
        """Keeps only datapoints inside the relevant window,
        or outside, if outside=True is passed.

        If remove=False, returns but doesn't remove the spikes.

        Returns: the indices of spikes and of clusters in the area."""
        (xmin, xmax, ymin, ymax) = rectangle

        dx, dy = self.__data
        numclus = self.NClusters()
        initialdata = self.NData()
        if not outside:
            condition = [x & y & z & w for (x, y, z, w) in zip(dx <= xmax,
                         dx >= xmin, dy <= ymax, dy >= ymin)]
        else:
            condition = [-(x & y & z & w) for (x, y, z, w) in zip(dx <= xmax,
                         dx >= xmin, dy <= ymax, dy >= ymin)]
        d_ind_kept = np.where(condition)[0]
        cid_kept_all = self.__ClusterID[d_ind_kept]
        c_ind_kept = np.unique(cid_kept_all)
        if remove:
            if numclus:
                newID = -np.ones(numclus, dtype=np.int)
                newID[c_ind_kept] = np.array(range(len(c_ind_kept)))
                # update temporarily the ClusterID vector
                self.__ClusterID = newID[self.__ClusterID]
                # delete data whose cluster was deleted, and clusters
                self.__ClusterID = self.__ClusterID[d_ind_kept]
                self.__c = self.__c[:, c_ind_kept]
            self.Backup()
            self.KeepOnly(d_ind_kept)
            print('Crop removed ' + str(numclus-self.NClusters()) +
                  ' clusters and ' + str(initialdata-self.NData()) +
                  ' datapoints.')
        return d_ind_kept, c_ind_kept

# UTILITY

    def UpdateExperimentIndices(self, myInds):
        """This is used when applying a filter, to keep track
        of the indices at which new stimulation protocols begin"""
        if len(self.__expinds) > 1:
            for n, i in enumerate(self.__expinds[1:]):
                self.__expinds[n + 1] = np.where(myInds >= i)[0][0]
            print('New experiment indices: '+str(self.__expinds))

    def KeepOnly(self, ind_kept):
        """This is used to remove datapoints that were filtered out
        and update the arrays. When the data are clustered, more
        updates need to be done"""
        # does not act on clusters!
        self.__data = self.__data[:, ind_kept]
        if np.size(self.__shapes):
            self.__shapes = self.__shapes[:, ind_kept]
        if np.size(self.__times):
            self.__times = self.__times[ind_kept]
        self.UpdateExperimentIndices(ind_kept)

    def Backup(self):
        """Creates a checkpoint, to be used for a subsequent
        call to UndoLast()"""
        self.__backup = {0: self.__data, 1: self.__ClusterID,
                         2: self.__c, 3: self.__shapes, 4: self.__times}

    def UndoLast(self):
        """The object restores the data as it was before
        the last call of a filter, or Backup()."""
        self.__data, self.__ClusterID, self.__c, self.__shapes, \
            self.__times = self.__backup[0], self.__backup[1], \
            self.__backup[2], self.__backup[3], self.__backup[4]

# OTHER
    def ClusterWidth(self):
        """Compute the spread of spikes within each cluster. This can be
        used to distinguish well centred clusters from those where spike
        locations were poorly mapped.

        Note this method uses the cluster centre locations returned by
        the mean shift algorithm, which differ from the centre of mass
        if PCA is used alonside position.

        Returns:
        clwidth -- The standard deviation of the spike locations relative
        to the cluster centre.
        """

        clwidth = np.zeros(self.NClusters())
        for n in range(self.NClusters()):
            centre = self.ClusterLoc()[:2, n]
            inds = self.__ClusterID == n
            clwidth[n] = np.std(np.sqrt((self.__data[0, inds]-centre[0])**2 +
                                        (self.__data[1, inds]-centre[1])**2))
        return clwidth

    def QualityMeasures(self, scorePCA=None, ncomp=None):
        return QualityMeasures(self, scorePCA, ncomp)

    def ShapeClassifier(self):
        return ShapeClassifier(self)


# A separate class to build a classifier.
class ShapeClassifier(object):
    def __init__(self, spikeobj):
        self.spikeobj = spikeobj

    def BadShapesByDensity(self, nbins=[64, 64],
                           percentile=0.5,
                           maxn=None,
                           min_thr=5,
                           normalise=False):
        l = self.spikeobj.Locations() + 0.5  # remove +0.5 at some point
        hg, bx, by = np.histogram2d(l[0], l[1], nbins)
        mindensity = np.min(hg[hg > 0])
        density_thr = np.max((np.percentile(hg.flatten(), percentile),
                             mindensity + min_thr))  # +5 is also arbitrary!
        binspanx = (np.max(l[0]) - np.min(l[0]))/nbins[0]*1.001
        binspany = (np.max(l[1]) - np.min(l[1]))/nbins[1]*1.001
        nbx = ((l[0] - np.min(l[0])) // binspanx).astype(int)
        nby = ((l[1] - np.min(l[1])) // binspany).astype(int)
        indbad = np.where(hg[nbx, nby] <= density_thr)[0]
        if maxn is not None:
            indbad = np.sort(np.random.permutation(indbad)[:maxn])
        if normalise:
            normed = lambda X: X/np.max(np.abs(X), axis=0)
            badshape = np.median(normed(self.spikeobj.Shapes()[:, indbad]),
                                 axis=1)
        else:
            badshape = np.median(self.spikeobj.Shapes()[:, indbad], axis=1)
        print("Working with " + str(len(indbad)) +
              " examples of bad shapes.")
        return badshape, indbad

    def GoodShapesByAmplitude(self, amp_thr, maxn=None, normalise=False):
        fakeampl = -np.min(self.spikeobj.Shapes(), axis=0)
        indgood = np.where(fakeampl > amp_thr)[0]
        if maxn is not None:
            indgood = np.sort(np.random.permutation(indgood)[:maxn])
        print("Working with " + str(len(indgood)) +
              " examples of good shapes.")
        if normalise:
            normed = lambda X: X/np.max(np.abs(X), axis=0)
            goodshape = np.median(normed(self.spikeobj.Shapes()[:, indgood]),
                                  axis=1)
        else:
            goodshape = np.median((self.spikeobj.Shapes()[:, indgood]), axis=1)
        return goodshape, indgood

    def FitClassifier(self, pcascores, indgood, indbad):
        # create a matrix of waveform PC projections
        pcs = np.hstack((pcascores[:, indbad], pcascores[:, indgood]))
        # the training labels
        # WHY do we use 0 and 1 instead of projections of some kind?
        labels = np.append(np.zeros(len(indbad)), np.ones(len(indgood)))
        # fit the classifier
        classifier = svm.SVC(kernel='rbf', class_weight='balanced')
        # use this for sklearn <0.17
        # classifier = svm.SVC(kernel='rbf',class_weight='auto')
        classifier.fit(pcs.T, labels)
        # get the labels for the whole data set
        score = classifier.predict(pcascores.T).astype(int)
        print("Classified as bad: "+str(np.sum(score == 0)) +
              ", and as good: "+str(np.sum(score == 1)))
        return score


class QualityMeasures(object):
    def __init__(self, spikeobj, scorePCA=None, ncomp=None):
        if np.size(spikeobj.ClusterID()) == 0:
            raise ValueError('No clustering was performed')
        self.spikeobj = spikeobj
        if scorePCA is None:
            scorePCA = self.spikeobj.ShapePCA(ncomp=ncomp, white=True)
        self.scorePCA = scorePCA

    def Neighbours(self, cl_idx, d, min_neigh_size=0):
        clocs = self.spikeobj.ClusterLoc()
        clsizes = self.spikeobj.ClusterSizes()
        dists = (clocs[0] - clocs[0, cl_idx])**2 + \
                (clocs[1] - clocs[1, cl_idx])**2
        return np.where((dists > 0) & (dists < d**2) &
                        (clsizes >= min_neigh_size))[0]

    def GaussianOverlapGroup(self, clnumbers, mode="both"):
        fourvector = np.vstack((self.spikeobj.Locations(),
                               self.scorePCA[:4, :]))
        fstd = np.std(fourvector, axis=1)
        fstd[:2] = 1
        fmean = np.mean(fourvector, axis=1)
        fmean[:2] = 0
        spLabels = self.spikeobj.ClusterID()
        inds = []
        for j in clnumbers:
            inds.append(np.where(spLabels == j)[0])
        data = []
        if mode == "both":
            for ind in inds:
                data.append((fourvector[:, ind].T-fmean)/fstd.T)
        elif mode == "XY":
            for ind in inds:
                data.append((fourvector[:2, ind].T-fmean[:2])/fstd[:2].T)
        elif mode == "PCA":
            for ind in inds:
                data.append((fourvector[2:, ind].T-fmean[2:])/fstd[2:].T)
        else:
            raise ValueError("Acceptable modes are 'all', 'PCA' and 'XY'")
        return self._gaussian_overlap(data)  # confusion matrix

    def _data_gaussian_overlap(p):
        '''
        Fit a len(p)-component Gaussian mixture model to a set of clusters,
        estimate the cluster overlap and return a matrix containing estimated
        false positives and negatives.

        Data is provided as list in p, each an array containing PCA projections
        or locations or both.

        This method is based on:
        Hill, Daniel N., Samar B. Mehta, and David Kleinfeld.
        Quality metrics to accompany spike sorting of extracellular signals.
        Journal of Neuroscience 31.24 (2011): 8699-8705.

        From the original description by Hill et al.:
        The percent of false positive and false negative errors are estimated
        for both classes and stored as a confusion matrix. Error rates are
        calculated by integrating the posterior probability of a
        misclassification. The integral is then normalized by the number of
        events in the cluster of interest.

        Returns:
        confusion - a confusion matrix, diagonals have fase positive, and
        off-diagonals false negatives
        '''

        ncl = len(p)
        # pre-compute means and covariance matrices
        estCent = np.array([np.median(p[i], axis=0) for i in range(len(p))])
        estCov = np.array([np.cov(p[i].T) for i in range(len(p))])
        g = mixture.GMM(n_components=ncl, covariance_type='full',
                        params='wmc', init_params='w',
                        min_covar=1e-2, tol=1e-3)
        g.means_ = np.vstack(estCent)
        g.covars_ = estCov
        data = np.concatenate(p)
        g.fit(data)
        if g.converged_ is False:
            print("not converged")

        # get posterior probabilities
        pr = []
        for i in range(ncl):
            pr.append(g.predict_proba(p[i]))

        # get indices in case clusters are mixed up
        # assign each GMM cluster to the nearest cluster the first two dims
        pInds = np.zeros(ncl, dtype=int)
        d = euclidean_distances(np.vstack(estCent)[:, :2], g.means_[:, :2])
        for i in range(ncl):
            ind = np.argmin(d)
            pInds[np.floor(ind/ncl)] = ind % ncl
            d[:, ind % ncl] = 10
            d[np.floor(ind/ncl)] = 10

        # compute the confusion matrix entries
        confusion = np.zeros((ncl, ncl))
        for i in range(ncl):
            # FP
            confusion[pInds[i], pInds[i]] = np.sum(
                np.mean(pr[i][:, np.setxor1d(i, range(ncl))], axis=0))
            # FNs
            for j in np.setxor1d(i, range(0, ncl)):
                confusion[pInds[i], pInds[j]] = np.sum(pr[j][:, i])/len(p[i])

        # if plotResult is True:  # plot the result
        #     ax = plt.gca()
        #     nShow = 500
        #     for i in range(ncl):
        #         plt.scatter(p[i][:nShow,pltDims[0]], p[i][:nShow,pltDims[1]],
        #        s=2,color=colors[2*i],label=str(np.round(confusion[i,i],3)))
        #     leg = plt.legend(loc=1, ncol=1,frameon=False,fontsize=14,
        #                      handlelength=0,handleheight=0,columnspacing=0.2,handletextpad=0)
        #     for i, txt in enumerate( leg.get_texts()):
        #         txt.set_color(colors[2*i])
        #     for item in leg.legendHandles:
        #         item.set_visible(False)
        #     plt.axis('equal')
        #     for i in range(ncl):
        #         v, w = np.linalg.eigh(g._get_covars()[pInds[i]][slice(
        #            pltDims[0],pltDims[1]+1),slice(pltDims[0],pltDims[1]+1)])
        #         u = w[0] / np.linalg.norm(w[0])
        #         angle = np.arctan2(u[1], u[0])
        #         angle = 180 * angle / np.pi  # convert to degrees
        #         v *= 9./2
        #         ell = mpl.patches.Ellipse(g.means_[pInds[i], pltDims], v[0],
        #            v[1], 180 + angle, color=colors[2*i])
        #         ell.set_clip_box(ax.bbox)
        #         ell.set_linewidth(0.7)
        #         ell.set_edgecolor('k')
        #         ell.set_alpha(0.55)
        #         ax.add_artist(ell)

        return confusion
