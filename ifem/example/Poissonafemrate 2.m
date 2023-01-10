%% RATE OF CONVERGENCE OF AFEM: P1 Linear Element for Poisson 
%
% This example is to show the rate of convergence of the lowest order
% finite element approximation of the second order elliptic equation.
%
% # Lshape problem.
% # Kellogg problem.

%% Lshape problem
[node,elem] = squaremesh([-1,1,-1,1],1);
[node,elem] = delmesh(node,elem,'x>0 & y<0');
bdFlag = setboundary(node,elem,'Dirichlet');
pde = Lshapedata;
format shorte
option.L0 = 1;
option.maxIt = 25;
option.printlevel = 1;
option.plotflag = 1;
err = afemPoisson(node,elem,pde,bdFlag,option);
figure;
showrate2(err.N,err.H1,10,'k-*','||Du-Du_h||',err.N,err.eta,10,'-+','eta');

%% Kellogg problem
[node,elem] = squaremesh([-1 1 -1 1], 0.5);
bdFlag = setboundary(node,elem,'Dirichlet');
pde = Kelloggdata;
option.L0 = 1;
option.maxIt = 100;
option.maxN = 1e4;
option.theta = 0.2;
option.plotflag = 1;
err = afemPoisson(node,elem,pde,bdFlag,option);
figure;
showrate2(err.N,err.H1,20,'k-*','||Du-Du_h||',err.N,err.eta,40,'-+','eta');
% latexerrtable(err.N,[err.H1 err.eta])