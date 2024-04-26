# AIUQ
This package is for scattering analysis of microscopy. The code have been tested in MATLAB version 2022b.

Description:

The package allows users to simulate 2D movements of some particels from different anisotropic or isotropic stochastic processes and model the mean squared displacement (MSD) given a dataset using Gaussian processes. Parameter estimation is performed in a fast approach, making use of the Toeplitz structure of covariance matrices.  

References:
  1. M. Gu, Y. He, X. Liu, and Y. Luo, Ab initio uncertainty quantification in scattering analysis of microscopy, arXiv preprint arXiv:2309.02468
  2. Y. Ling, Superfast Inference for Stationary Gaussian Processes in Particle Tracking Microrheology, Ph.D. thesis, University of Waterloo (2019).
  3. Cerbino, R. and Trappe, V., 2008. Differential dynamic microscopy: probing wave vector dependent dynamics with a microscope. Physical review letters, 100(18), p.188102.

Installation:

To use this package, please install the optimization toolbox in MATLAB and install the FFTW3 library. For Windows users and Mac users whose system is based on  X86_64 architecture (usually before 2021), FFTW3 library can be linked easily, and for Mac users with ARM architecture (M chips) the procedure may be more complicated.

To install and link FFTW3 to your system, please refer to the FFTW_install_tutorial_2024.pdf


Contents:

example_anisotropic.m - Show anisotropic examples of some simulations and implement aniso_SAM() to real intensity data
aniso_simulation.m - Simulation code for different anisotropic stochastic processes
anisotropic_neg_log_lik - Compute the negative log-likelihood of the processed intensity data. 
anisotropic_log_lik_grad - Compute the derivative of negative log-likelihood w.r.t parameters. 
aniso_SAM - Parameter estimation 
example.m - Show isotropic examples of some simulations and implement aniso_SAM() to real intensity data
simulation.m - Simulation code for different isotropic stochastic processes
neg_log_lik - Compute the negative log-likelihood of the processed intensity data. 
log_lik_grad - Compute the derivative of negative log-likelihood w.r.t parameters. 
SAM - Parameter estimation 
cpp files: log_density, phi_log_det, toeplitz_prod, toeplitz_solve, toeplitz_trace_gradient, toeplitz_trace_hess

The main functions are aniso_simulation(), simulation() aniso_SAM() and SAM().

Authors:
Xubo Liu, Yue He, Mengyang Gu

Department of Statistics and Applied Probability, University of California, Santa Barbara

Email: mengyang@pstat.ucsb.edu, xubo@ucsb.edu, yuehe@ucsb.edu
