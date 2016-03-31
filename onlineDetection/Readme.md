# Online spike detection

Developed and written by [Oliver Muthmann](ollimuh@googlemail.com).

Reference: J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.

Note this implementation is currently only compatible with the old .brw file format used by 3Brain to store Biocam4096 recordings with Brainwave 2.0. A generic version with a python interface will be ready shortly.

## Installation and use

1. Download the current *BrwExtReader.dll* from [http://www.3brain.com/index.php/5/Downloads](http://www.3brain.com/index.php/5/Downloads) or use the version provided here
2. Compile Main.cs (use any C# compiler for that, e.g. Monodevelop on Linux), or use the executable provided
3. Run (e.g. using Mono)

## Implementation

The online detection will produce a single text file containing (integers arranged in three columns) channel ids, time stamps (in units of frames) and amplitudes (in units of v/64).

Some parameters can be adjusted by the experimenter. The repolarisation threshold ensures that the measured voltage goes back towards baseline. At least one frame during the maximum width of the spike should be larger than that value. Note, however, that this can be strongly distorted by noise and therefore one should not use a high value here.

The output from the spike detection is then read and processed by the PYTHON tools provided in [postProcessing](../postProcessing).
