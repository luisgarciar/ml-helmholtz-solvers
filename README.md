# Multilevel solvers for the Helmholtz equation

This repository contains an implementation of the multilevel methods for the solution of linear systems resulting from the discretization of Helmholtz problems described in the paper [GRN20] and the PhD dissertation [GR22]. It is based on a custom version of the [iFEM](https://github.com/lyc102/ifem) library [C09], and also contains standalone implementations of the multigrid method and the shifted Laplace preconditioner [EOV06] for Helmholtz problems discretized with finite differences and finite elements. In addition to solving the problems presented in the paper and the dissertation, the code in the repository can be adapted to solve other Helmholtz problems in 2D.  


## Installation and Usage
The software has been developed in MATLAB 2017b and tested with versions up to MATLAB 2020a. To clone this repository, navigate using the terminal to your desired location and type
`git clone https://github.com/luisgarciar/ml_helmholtz_solvers.git`

Next, add the path to `ml_helmholtz_solvers` to your MATLAB path, which can be done using the graphical interface in MATLAB (File -> Set Path -> Add with Subfolders).

## Examples and Tests
The folders `FD` and `FEM` contain the implementations of the solver for finite difference and finite element discretizations respectively. Each of these folders contains tests and examples showcasing the various functionalities of the solvers. The folder  `num_exp` contains a variety of numerical experiments for Helmholtz problems in 1D and 2D. The experiments in the paper [GR22] are located in the folder `exp_TL_paper` .


## References

[GR22] L. García Ramos, [Polynomial and multilevel preconditioners for the Helmholtz equation based on the shifted Laplacian](https://depositonce.tu-berlin.de/handle/11303/17048). PhD Thesis, Technische Universität Berlin, 2022.

[GRN20 ] L. García Ramos and R. Nabben, [A two-level shifted Laplace preconditioner for Helmholtz problems: Field-of-values analysis and wavenumber-independent convergence](https://arxiv.org/abs/2006.08750). arXiv: 2006.08750.

[EOV06] Y. A. Erlangga, CW Oosterlee, C Vuik
[A novel multigrid based preconditioner for heterogeneous Helmholtz problems](https://epubs.siam.org/doi/10.1137/040615195). SIAM Journal on Scientific Computing 27 (4), 1471-1492.

[C09] L. Chen, [iFEM: an integrated finite element method package in MATLAB](https://www.math.uci.edu/~chenlong/Papers/iFEMpaper.pdf). Technical Report, University of California at Irvine, 2009.
