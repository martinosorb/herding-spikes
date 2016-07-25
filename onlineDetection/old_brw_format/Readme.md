# Online spike detection

Developed and written by [Oliver Muthmann](ollimuh@googlemail.com). Cython port by [Albert Puente Encinas](https://github.com/albertpuente) and [Matthias Hennig](http://homepages.inf.ed.ac.uk/mhennig/index.html).

Reference: J.-O. Muthmann, H. Amin, E. Sernagor, A. Maccione, D. Panas, L. Berdondini, U.S. Bhalla, M.H. Hennig MH (2015). [Spike detection for large neural populations using high density multielectrode arrays](http://journal.frontiersin.org/article/10.3389/fninf.2015.00028/abstract). Front. Neuroinform. 9:28. doi: 10.3389/fninf.2015.00028.

This branch provides a version which can process raw data from any source.

This original implementation is compatible with the old .brw file format used by 3Brain to store Biocam4096 recordings with Brainwave 2.0.

## Installation and use

1. Download the current *BrwExtReader.dll* from [http://www.3brain.com/index.php/5/Downloads](http://www.3brain.com/index.php/5/Downloads) or use the version provided here
2. Use the executable provided, or compile Main.cs (use any C# compiler for that, e.g. [gmcs](http://www.mono-project.com/docs/about-mono/languages/csharp/) or [Monodevelop](http://www.monodevelop.com/) on Linux).
3. Run (e.g. using [Mono](http://www.mono-project.com/)), a graphical user interface will appear.
