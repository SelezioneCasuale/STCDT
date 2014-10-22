STCDT
=====

## The Spike Train Community Detection Toolbox

This repo hosts the current usable version of the Spike Train Community Detection toolbox developed by the Mark Humphries lab at the University of Manchester. Code included here implements the algorithms from _Humphries, M. D. (2011). Spike-train communities: finding groups of similar spike trains. The Journal of Neuroscience, 31(6), 2321-2336_ . A pdf of the paper is provided in this repo for reference.

Matlab, Python, and Spark (notebook format, very developmental) versions of the code are provided.

Each version of the code is at a different stage of development.

### Matlab:
Most stable. Version 1.2 released on 12.6.2011. Written by Mark Humphries. The code here mirrors that available [here](http://www.systemsneurophysiologylab.ls.manchester.ac.uk/code/analysis/ "Humphries lab website").

### Python:
A port of the matlab code to Python. It has been tested on the synthetic data given in the matlab code, but hasn't been used "in anger" on larger datasets.

### Spark:
IPython notebook for a large scale implementation of the algorithms, currently in development by Mathew Evans. The code was built on top of a python script written by Thomas Sharp @ Riken institute. The notebook can be viewed statically here : [http://nbviewer.ipython.org/github/mathewzilla/STCDT/blob/master/Spark/Community_Detection.ipynb](http://nbviewer.ipython.org/github/mathewzilla/STCDT/blob/master/Spark/Community_Detection.ipynb)


#### Getting in touch
If you have any enquiries about any part of the code contact me (Mat Evans) either here through github or via email [mathew.evans@manchester.ac.uk](mailto:mathew.evans@manchester.ac.uk "I have the internets")
