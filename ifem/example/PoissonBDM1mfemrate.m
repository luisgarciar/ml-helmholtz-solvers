%% RATE OF CONVERGENCE OF MIXED FINITE ELEMENT METHOD (BDM1-P0) FOR POISSON EQUATION
%
% This example is to show the rate of convergence of mixed finite element
% (RT0-P0) approximation of the Poisson equation on the unit square:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^2$$
%
% for the following boundary conditions:
%
% # Pure Dirichlet boundary condition. $\Gamma _D = \partial \Omega$. 
% # Pure Neumann boundary condition. $\Gamma _N = \partial \Omega$.
% # Mix Dirichlet and Neumann boundary condition. $u=g_D \hbox{ on }\Gamma_D, \quad \nabla u\cdot n=g_N \hbox{ on }\Gamma_N.$
%
% Written by Ming Wang and improved by Long Chen.

%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
pde = sincosNeumanndata;
option.L0 = 2;
option.maxIt = 4;
option.printlevel = 1;
option.elemType = 'BDM1';

%% Non-empty Dirichlet boundary condition.
option.plotflag = 1;
bdFlag = setboundary(node,elem,'Dirichlet');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Pure Neumann boundary condition.
option.solver = 'tripremixpoisson';
bdFlag = setboundary(node,elem,'Neumann');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Mix Dirichlet and Neumann boundary condition.
option.solver = 'uzawapcg';
bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
mfemPoisson(node,elem,pde,bdFlag,option);

%% Conclusion
%
% The optimal rates of convergence for u and sigma are observed, namely,
% 1st order for L2 norm of u, H(div) norm of sigma, and 2nd order for 
% L2 norm of sigma. 
%
% The 2nd order convergent rates between two discrete functions ||uI-uh|| 
% is known as superconvergence.
%
% The 3rd order convergent rates between two discrete functions 
% ||sigmaI-sigmah|| for the Pure Neumman boundary condition 
% is known as superconvergence. The 2nd accurate numerical quadrature
% is required for the integral of rhs to observe such superconvergence.
%
% Triangular preconditioned GMRES and Uzawa preconditioned CG converges
% uniformly in all cases. Traingular preconditioner is two times faster than
% PCG although GMRES is used.