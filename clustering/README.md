Spike sorting
=============

This code performs fully automated spike sorting on events detected and localised with the interpolating spike detection tools in  the sub-project [postProcessing](../postProcessing).


The full workflow is illustrated in a jupyter notebook: [Sorting_workflow_demonstration.ipynb](Sorting_workflow_demonstration.ipynb)

Sorting is very fast, these are run times measured on 10 cores Intel Xeon E5-2630, 2.60GHz:

Number of spikes | Run time (h:mm:ss)
------------------|---------
10000000 | 0:06:27
6000000 | 0:03:34
2000000 | 0:01:58

## Clustering very large data sets

The parallelised mean-shift implementation in sklearn fails when the number of seeds is very large (many tens of millions of spikes).

This problem is solved in an alternative implementation, where the seeds are broken down into batches, and these batches are executed in parallel. This may generally me more efficient as the overhead is reduced. To use this, replace mean_shift_.py with [this file](mean_shift_.py). The file is part of the sklearn distribution, if installed locally, it can be found in `~/.local/lib/python3.4/site-packages/sklearn/cluster/`.

## References

A description of the methods is in this paper:

G. Hilgen, M. Sorbaro, S. Pirmoradian, J.-O. Muthmann, I. Kepiro, S. Ullo, C. Juarez Ramirez, A. Puente Encinas, A. Maccione, L. Berdondini, V. Murino, D. Sona, F. Cella Zanacchi, E. Sernagor, M.H. Hennig (2016). [Unsupervised spike sorting for large scale, high density multielectrode arrays.](http://www.cell.com/cell-reports/fulltext/S2211-1247(17)30236-X) Cell Reports 18, 2521–2532. bioRxiv doi: [http://dx.doi.org/10.1101/048645](http://dx.doi.org/10.1101/048645).

This method was developed by [Martino Sorbaro](http://martinosorb.github.io) and [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html).
