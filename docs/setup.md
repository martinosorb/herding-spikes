How to get started
------------------

## Requirements

Make sure the following is on your system:

* Python packages:
  * [h5py](http://www.h5py.org/)
  * [numpy](http://www.numpy.org/)
  * [matplotlib](http://matplotlib.org/)
  * [sklearn](http://scikit-learn.org/stable/)
  * [cython](http://cython.org/)


* Compilers
  * Linux: C++11 compiler (gcc/g++)
  * Windows: [Microsoft Visual C++ Compiler for Python 2.7 ](https://www.microsoft.com/en-us/download/details.aspx?id=44266) or any full-fledged Microsoft Visual C++ compiler (>= 2008).

In the following it assumed a linux shell (or similar) is used, some operations will of course differ under Windows.

## Getting the code

```
git clone -b detection-cpp https://github.com/martinosorb/herding-spikes.git
```

Note here we clone the branch ```detection-cpp```, which has the (beta, but very fast and hdf5 compatible) C++/Cython based detection code.

## Compile the detection code

Two different spike detection methods are available, both implemented in C++ and with a nice Cython interface. To compile these, use:

```
~/herding-spikes> cd herding-spikes/onlineDetection
~/herding-spikes/onlineDetection> python setup.py build_ext --inplace
~/herding-spikes/onlineDetection> cd ../interpolatingDetection/
~/herding-spikes/interpolatingDetection> python setup.py build_ext --inplace
```

All code is (should be) compatible with both python version 2 and 3.

## Analysing data

Move on to the [spike sorting tutorial](sorting-tutorial.md).
