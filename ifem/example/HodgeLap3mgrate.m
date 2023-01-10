
%%
close all;  clear all
option.elemType = 'ND0';

%% Options
option.maxIt = 4;
% option.mg.solver = 'VCYCLE';    % V,W;
% option.mg.coarsematrix = 'G'; % Galerkin formulation
option.mg.smoothingstep = 3;    % Smoothing step.
option.mg.smoothingratio = 1.5; % ratio of variable smoothing
option.mg.Vit = 1;      % number of cycles for Schur complement   

%% Cube
pde = HodgeLaplacian3Edata1;
[node,elem] = cubemesh([0,1,0,1,0,1],0.25);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% Diagonal Preconditioner
option.solver = 'diag';
mfemHodgeLap3(node,elem,pde,bdFlag,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap3(node,elem,pde,bdFlag,option);

%% Lshape domain
pde = fveconedata;
[node,elem] = cubemesh([-1,1,-1,1,-1,1],0.5);
% [node,elem] = delmesh(node,elem,'x>0 & y<0');
[node,elem] = delmesh(node,elem,'x<0 & y<0 & z>0');
showboundary3(node,elem);
bdFlag = setboundary3(node,elem,'Dirichlet');

%% Diagonal Preconditioner
option.solver = 'diag';
mfemHodgeLap3(node,elem,pde,bdFlag,option);

%% Triangular Preconditioner
option.solver = 'tri';
mfemHodgeLap3(node,elem,pde,bdFlag,option);