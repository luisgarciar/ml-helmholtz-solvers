%% RATE OF CONVERGENCE OF CUBIC FINITE ELEMENT METHOD
%
% This example is to show the rate of convergence of linear finite element
% approximation of the Poisson equation on the unit square:
%
% $$- \Delta u = f \; \hbox{in } (0,1)^2$$
%
% for the following boundary condition:
%

%% Set up problem
% PDE and Boundary condition.
pde = simpledata; % f = 1, g_D = 0
% FEM
option.elemType = 'P3';

%% Options
option.L0 = 0;
option.maxIt = 4;
option.printlevel = 1;
option.plotflag = 1;

%% Case 1: Triangle
[node,elem] = regpolygon(3,0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
showmesh(node,elem);
femPoisson(node,elem,pde,bdFlag,option);

%% Case 2: Square
[node,elem] = squaremesh([0,1,0,1],0.25); 
bdFlag = setboundary(node,elem,'Dirichlet');
showmesh(node,elem);
femPoisson(node,elem,pde,bdFlag,option);

%% Case 3: Pentagon
[node,elem] = regpolygon(5,0.5);
showmesh(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
femPoisson(node,elem,pde,bdFlag,option);

%% Case 4: Hexagon
[node,elem] = regpolygon(6,0.5);
showmesh(node,elem);
bdFlag = setboundary(node,elem,'Dirichlet');
femPoisson(node,elem,pde,bdFlag,option);