%Test of functions mg_setup, lin_interpol, fwrestriction
clearvars -global
close all;
clc;

k    = 40;          %wavenumber
ppw  = 20;          %min points per wavelength%
npcc  = 3;          %number of points in coarsest grid
bc   = 'som';       %boundary conditions
dim  =  2;          %dimension
eps  = 0.5*k^2 ;    %imaginary shift of shifted Laplacian

%number of fine points and levels
%computed from number of points in 
%coarsest grid and wavenumber
[npf,numlev] = fd_npc_to_npf(npcc,k,ppw);  

% Test of V-cycle on a Helmholtz or Poisson problem
% Construct matrices on all grids and interpolation operators
op_type = 'rd';

[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npcc,numlev,bc,dim);
A = mg_mat{1}; N = length(A);
b  = zeros(N,1);
x0 = randn(N,1);

% Parameters of V-cycle and Jacobi iteration
npre = 1; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 10;
r0   = norm(b-A*x0);

%% Test of multigrid on Helmholtz problem
%Use the profiler to evaluate performance of code
profile on
[x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
profile off

r0 
r1 = norm(b-A*x_sol)/norm(b-A*x0) %Multigrid diverges on Helmholtz problems!


%% Test of multigrid on shifted Laplacian problem
% [mg_mat,mg_split,restrict,interp] = mg_setup(M,k,eps,op_type,numlev,bc,dim);
% 
% profile on
% [x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
% profile off
% 
% r0
% r1 = norm(b-M*x_sol)/norm(b-M*x0)

