# BrainPalimpsest

_What else than a natural and mighty palimpsest is the human brain?_ 
_-- Thomas De Quincey_

This toolbox provides deconvolution methods for BOLD-fMRI data and yields simultanesouly  
spatiotemporal images of:

- underlying neural activity,
- neuroglial drive, 
- cerebral blood flow (CBF), 
- cerebral blood volume (CBV),
- deoxygenated hemoglobin (dHb) concentration, and, 
- response modes  W (wave), L (local-oscillating), and D (local-decaying).


Some important aspects before using the toolbox:

1. Install Freesurfer from https://surfer.nmr.mgh.harvard.edu/

2. [Build a `mex` file](https://au.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html) from `inplaceprod.c `. This step requires to have [a C Matlab-compatible compiler installed](https://au.mathworks.com/support/compilers.html). 
