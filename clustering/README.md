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

A description of the methods is availbe in this preprint:

G. Hilgen, M. Sorbaro, S. Pirmoradian, J.-O. Muthmann, I. Kepiro, S. Ullo, C. Juarez Ramirez, A. Maccione, L. Berdondini, V. Murino, D. Sona, F. Cella Zanacchi, U. Bhalla, E. Sernagor, M.H. Hennig (2016). [Unsupervised spike sorting for large scale, high density multielectrode arrays.](http://dx.doi.org/10.1101/048645) bioRxiv doi: http://dx.doi.org/10.1101/048645.

This method was developed by [Martino Sorbaro](http://martinosorb.github.io) and [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html).
