%% RATE OF CONVERGENCE OF MIXED FINITE ELEMENT METHOD (RT0-P0) FOR DARCY'S EQUATIONS
%
% This example is to show the rate of convergence of mixed finite element
% (RT0-P0) approximation of the Darcy's equations.


close all
%% Setting
[node,elem] = squaremesh([0,1,0,1],0.25); 
[bdNode,bdEdge,isBdNode] = findboundary(elem);
isIntNode = ~isBdNode;
node(isIntNode,:) = node(isIntNode,:) + 0.25*0.5*rand(sum(isIntNode),2);
% option.L0 = 1;
option.maxIt = 5;
option.printlevel = 1;
option.elemType = 'RT0';
option.refType = 'red';

%% Jump tensor
pde = Darcydata3;
% perturbe the grid
showmesh(node,elem);
% option.solver = 'uzawapcg';
% option.solver = 'tripremixpoisson';
option.solver = 'mg';
% bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
% bdFlag = setboundary(node,elem,'Dirichlet');
bdFlag = setboundary(node,elem,'Neumann');
mfemDarcy(node,elem,pde,bdFlag,option);

%% Isotropic tensor
% pde = Darcydata1;
pde = Darcydata2;
% option.solver = 'uzawapcg';
% option.solver = 'tripremixpoisson';
option.solver = 'mg';
% bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
% bdFlag = setboundary(node,elem,'Dirichlet');
bdFlag = setboundary(node,elem,'Neumann');
mfemDarcy(node,elem,pde,bdFlag,option);

%% Anisotropic tensor
pde = Darcydata2;
% option.solver = 'uzawapcg';
option.solver = 'tripremixpoisson';
% bdFlag = setboundary(node,elem,'Dirichlet','~(x==0)','Neumann','x==0');
bdFlag = setboundary(node,elem,'Dirichlet');
mfemDarcy(node,elem,pde,bdFlag,option);
