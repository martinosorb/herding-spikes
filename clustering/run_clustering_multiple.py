from herdingspikes import *
import os
import sys
import glob
import getopt
import h5py
from datetime import datetime
from sklearn import svm

if __name__ == "__main__":

    def printargs():
        print('run_clustering_multiple.py -i <inputfile or file mask> [-m <mbf> -a <0/1> -p <0/1>]')

    if len(sys.argv)==1:
        printargs()
        sys.exit(2)

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:m:a:p:", [
                                   "ifile=", "mbf=", "align=", "presort="])
    except getopt.GetoptError:
        printargs()
        sys.exit(2)

    # defaults
    mbf = 20
    # set to True is spike shapes should be re-aligned
    align = False
    # set to True if spikes should be pre-filtered, advisable at low samplig
    # rates
    presortSpikes = True
    # set to True if shapes should not be saved (smaller file)
    noShapes = False
    # relative weight of PCA in clustering
    alpha = 0.32
    # kernel size
    h = 0.3
    # for presorting, spikes with larger amplitude will be used to train the classifier
    # for true spikes
    ampthreshold = 5
    # shapes offsets...
    shoffset = 0  # 14 for 21kHz 18 for 24khz 0 for 7khz
    upto = 0  # 50 for 21kHz 0 for 7khz

    for opt, arg in opts:
        print(opt, arg)
        if opt == '-h':
            printargs()
            sys.exit()
        elif opt in ("-i", "--ifile"):
            filemask = arg
        elif opt in ("-m", "--mbf"):
            mbf = int(arg)
        elif opt in ("-a", "--align"):
            align = bool(int(arg))
        elif opt in ("-p", "--presort"):
            presort = bool(int(arg))

    try:
        sfs = glob.glob(filemask)
        print('Sorting file(s):')
        for f in sfs:
            print(f)
    except:
        print('cannot find ' + filemask)
        sys.exit(2)

    try:
        f = h5py.File(sfs[0], 'r')
        samplig = f['Sampling'].value
        f.close()
        print('Sampling rate:' + str(samplig))
        if samplig < 8000: # set some defaults here
            shapeLenth = 22
        else:
            shapeLenth = 67
        print('Setting shape length to '+str(shapeLenth)+', check if this is correct!')
    except:
        print('error reading ' + sfs[0])

    if(len(sfs) > 1):
        multiFile = True
    else:
        multiFile = False

    # check file type and load the data (clumsy)
    f = h5py.File(sfs[0],'r')
    keys = f.keys()
    if not any([k=='AmplitudeThresholds' for k in keys]):
        print('Loading from previously clustered file(s)')
        if multiFile == True:
            O = LoadMultipleClustered(sfs, shapesrange=[0,26])
        else:
            O = spikeclass(sfs[0])
    else:
        if(len(sfs) > 1):
            O = ImportInterpolatedList(sfs)
        else:
            O = ImportInterpolated(sfs[0])
    f.close()

    outpostfix = '_clusteredX_' + str(h) + '_' + str(alpha) + '_' + str(mbf)
    if align == True:
        outpostfix = outpostfix + '_align'
    if multiFile == True:
        outpostfix = outpostfix + '_multi'
    if presortSpikes == True:
        outpostfix = outpostfix + '_presorted2'
    if noShapes == True:
        outpostfix = outpostfix + '_noshapes'
    print('writing with postfix: ' + outpostfix)

    # align shapes
    if align:
        print("Aligning shapes...",)
        O.AlignShapes()
        print("done.")

    # presort, remove noise
    if presortSpikes == True:
        print("Selecting good spikes...",)
        shcl = ShapeClassifier(O)
        badshape, indbad = shcl.BadShapesByDensity(nbins=[64, 64], percentile=0.5, maxn=2000, min_thr=5, normalise=False)
        goodshape, indgood = shcl.GoodShapesByAmplitude(ampthreshold, maxn=2000, normalise=False)
        # print('Spikes in both (should be 0):', np.intersect1d(indbad, indgood))
        scorePCA = O.ShapePCA(ncomp=6, white=True, offset=shoffset, upto=upto)
        score = shcl.FitClassifier(scorePCA, indgood, indbad)

        # plot results
        l = O.Locations()
        fig = plt.figure(figsize=(10, 4 * 3))
        nbins = [64 * 2, 64 * 2]
        cmap = plt.cm.hot
        cmap.set_bad('k')

        def formatplot(title):
            plt.colorbar()
            plt.axis('equal')
            plt.xlim((0, 65))
            plt.ylim((0, 65))
            ax.set_xlim(ax.get_xlim()[::-1])
            plt.title(title)

        ax = plt.subplot(321)
        hg, bx, by = np.histogram2d(l[0], l[1], nbins)
        rateMasked = np.ma.array(hg, mask=(hg == 0))
        with np.errstate(divide='ignore'):
            plt.pcolor(bx, by, np.log10(rateMasked), cmap=cmap, vmin=0, vmax=3)
        formatplot('all spikes')

        ax = plt.subplot(322)
        inds = np.where(score == 1)[0]
        hg1, bx, by = np.histogram2d(l[0][inds], l[1][inds], nbins)
        rateMasked = np.ma.array(hg1, mask=(hg1 == 0))
        with np.errstate(divide='ignore'):
            plt.pcolor(bx, by, np.log10(rateMasked), cmap=cmap, vmin=0, vmax=3)
        formatplot('good spikes')

        ax = plt.subplot(324)
        inds = np.where(score == 0)[0]
        hg12, bx, by = np.histogram2d(l[0][inds], l[1][inds], nbins)
        rateMasked = np.ma.array(hg12, mask=(hg12 == 0))
        with np.errstate(divide='ignore'):
            plt.pcolor(bx, by, np.log10(rateMasked), cmap=cmap, vmin=0, vmax=3)
        formatplot('bad spikes')

        ax = plt.subplot(326)
        rateMasked = np.ma.array(hg12 / hg, mask=(hg == 0))
        with np.errstate(divide='ignore'):
            plt.pcolor(bx, by, rateMasked, cmap=cmap, vmin=0, vmax=1)
        formatplot('fraction bad')

        ax = plt.subplot(325)
        plt.plot(badshape, label='bad')
        plt.plot(goodshape, label='good')
        plt.legend()

        ax = plt.subplot(323)
        nShow = 7000
        inds = np.where(score == 0)[0]
        plt.scatter(scorePCA[0, inds[:nShow]], scorePCA[
                    1, inds[:nShow]], c='r', s=6)
        inds = np.where(score == 1)[0]
        plt.scatter(scorePCA[0, inds[:nShow]], scorePCA[
                    1, inds[:nShow]], c='b', s=6)

        fname = sfs[0].replace('.hdf5', outpostfix) + '.png'
        plt.savefig(fname)
        plt.close()

        O.KeepOnly(np.where(score == 1)[0])
        print("done")

    # cluster data
    scorePCA = O.ShapePCA(ncomp=2, white=True, offset=shoffset, upto=upto)
    startTime = datetime.now()
    O.CombinedMeanShift(h, alpha, scorePCA, mbf=mbf)
    print('Time taken for sorting: ' + str(datetime.now() - startTime))

    # save data
    startTime = datetime.now()
    if multiFile == False:
        outfile = sfs[0].replace('.hdf5', outpostfix + '.hdf5')
        print('writing: '+outfile)
        if noShapes == False:
            O.Save(outfile, 'lzf')
        else:
            g = h5py.File(outfile, 'w')
            g.create_dataset("data", data=O.Locations())
            g.create_dataset("centres", data=O.ClusterLoc())
            g.create_dataset("cluster_id", data=O.ClusterID())
            g.create_dataset("times", data=O.Times())
            # g.create_dataset("shapes",data=O.Shapes()[:,inds],compression='lzf')
            g.create_dataset("Sampling", data=O.Sampling())
            g.close()
    else:
        for c, name in enumerate(sfs):
            ofname = name.replace('.hdf5', outpostfix + '.hdf5')
            print('writing: '+ofname)
            inds = O.ExperimentIndices(c)
            g = h5py.File(ofname, 'w')
            g.create_dataset("data", data=O.Locations()[:, inds])
            # g.create_dataset("expinds",data=self.__expinds)
            g.create_dataset("centres", data=O.ClusterLoc())
            g.create_dataset("cluster_id", data=O.ClusterID()[inds])
            g.create_dataset("times", data=O.Times()[inds])
            if noShapes == False:
                g.create_dataset("shapes", data=O.Shapes()[
                                 :, inds], compression='lzf')
            g.create_dataset("Sampling", data=O.Sampling())
            g.close()
    print('Time for saving: ' + str(datetime.now() - startTime))
