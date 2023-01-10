%% MULTIGRID OF FOR THE STOKES EQNS IN 2D
%
% This example is to show the convergence of multigrid methods for various
% finite element approximation of the Stokes equation on the unit square:
%
%       -div(mu*grad u) + grad p = f in \Omega,
%                        - div u = 0  in \Omega,
%                              u = g_D  on \Gamma.
%
% with the pure Dirichlet boundary condition.
%

%% Setting
clear all; close all;
% [node,elem] = squaremesh([0,1,0,1],0.25);
[node,elem] = circlemesh(0,0,1,0.25); option.mesh = 'circle';
% load Lshapemesh
% load Lshapeunstructure
% load flowpastcylindermesh
% node = [0 0; 1 0; 1 1];
% elem = [2 3 1];
% [node,elem] = uniformrefine(node,elem);
% [node,elem] = uniformrefine(node,elem);
showmesh(node,elem);
% pde = Stokesdata1; 
pde = StokesZulehnerdata;
bdFlag = setboundary(node,elem,'Dirichlet');
option.L0 = 1;
option.maxIt = 3;
option.printlevel = 1;
option.plotflag = 0;
option.rateflag = 0;

%% MG options
option.solver = 'mg';
% option.printlevel = 2;
option.smoothingstep = 2;
% option.smootherbarSp = 'SGS';

%% RT0-P0
% display('RT0-P0')
% option.elemType = 'RT0-P0';
option.elemType = 'BDM1B-P0';
% option.refType = 'bisect';
femStokesHdiv(node,elem,pde,bdFlag,option);

%% BDM1B-P0
% display('BDM1B-P0')
% option.elemType = 'BDM1B-P0';
% femStokesHdiv(node,elem,pde,bdFlag,option);