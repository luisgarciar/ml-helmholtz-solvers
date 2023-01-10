%% RATE OF CONVERGENCE OF CUBIC ELEMENT FOR POISSON EQUATION
%
% This example is to show the rate of convergence of cubic finite element
% approximation of the Poisson equation on the unit square:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^2$$
%
% for the following boundary condition:
%
% # Non-empty Dirichlet boundary condition. $u=g_D \hbox{ on }\Gamma_D, \quad \nabla u\cdot n=g_N \hbox{ on }\Gamma_N.$
% # Pure Neumann boundary condition. $\nabla u\cdot n=g_N \hbox{ on } \partial \Omega$.
% # Robin boundary condition. $g_R u + \nabla u\cdot n=g_N \hbox{ on }\partial \Omega$


%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
option.L0 = 1;
option.maxIt = 4;
option.maxN = 1e6;
option.printlevel = 1;
option.elemType = 'P3';

%% Non-empty Dirichlet boundary condition.
pde = sincosdata;
bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
femPoisson(node,elem,pde,bdFlag,option);

%% Pure Neumann boundary condition.
pde = sincosNeumanndata;
bdFlag = setboundary(node,elem,'Neumann');
femPoisson(node,elem,pde,bdFlag,option);

%% Pure Robin boundary condition.
option.plotflag = 0;
pdeRobin = sincosRobindata;
bdFlag = setboundary(node,elem,'Robin');
femPoisson(node,elem,pdeRobin,bdFlag,option);

%% Conclusion
%
% The optimal rate of convergence of the H1-norm (3rd order) and L2-norm
% (4th order) is observed. The order of ||DuI-Duh|| is 3rd order and
% thus no superconvergence exists between nodal interpolate and uh.
%
% MGCG converges uniformly in all cases.
