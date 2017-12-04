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

# Before using the toolbox

Some important aspects you need to do before using the toolbox:

1. Install Freesurfer from https://surfer.nmr.mgh.harvard.edu/

2. [Build a `mex` file](https://au.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html) from `inplaceprod.c `. This step requires to have [a C Matlab-compatible compiler installed](https://au.mathworks.com/support/compilers.html). 

# Citation

If you use our code in your research, please cite us as follows (will be updated as soon as the paper gets accepted):

J.C. Pang, K.M. Aquino, P.A. Robinson, T.C. Lacy, M.M. Schira, New fMRI windows reveal a palimpsest of brain activity, hemodynamics, and physiology, submitted to NeuroImage, 2017.

# Other references related to this work

J.C. Pang, P.A. Robinson, K.M. Aquino, N. Vasan, Effects of astrocytic dynamics on spatiotemporal hemodynamics: Modeling and enhanced data analysis, NeuroImage 147, 114-153, 2017.

T.C. Lacy, K.M. Aquino, P.A. Robinson, M.M. Schira, Shock-like haemodynamic responses induced in the primary visual cortex by moving visual stimuli, Journal of the Royal Society Interface, 13(125), 20160576, 2016.

J.C. Pang, P.A. Robinson, K.M. Aquino, Response-mode decomposition of spatio-temporal haemodynamics, Journal of the Royal Society Interface 13(118), 20160253, 2016.

K.M. Aquino, P.A. Robinson, M.M. Schira, M.J. Breakspear, Deconvolution of neural dynamics from fMRI data using a spatiotemporal hemodynamic response function, NeuroImage, 94, 203-215, 2014.

K.M. Aquino, P.A. Robinson, P.M. Drysdale, Spatiotemporal hemodynamic response functions derived from physiology, Journal of Theoretical Biology, 347(1), 118-136, 2014.

K.M. Aquino, M.M. Schira, P.A. Robinson, P.M. Drysdale, M.J. Breakspear, Hemodynamic traveling waves in human visual cortex, PLoS Computational Biology, 8(3), e1002435, 2012.

P.M. Drysdale, J.P. Huber, P.A. Robinson, K.M. Aquino, Spatiotemporal BOLD dynamics from a poroelastic hemodynamic model, Journal of Theoretical Biology, 265(4), 524-534, 2010.
