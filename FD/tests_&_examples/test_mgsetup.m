%Test of functions mg_setup, lin_interpol, fwrestriction
clear all; 
close all;
clc;

k    = 20;         %wavenumber
ppw  = 20;          %min points per wavelength%
npc  = 3;           %number of points in coarsest grid
bc   = 'dir';       %boundary conditions
dim  =  1;          %dimension
eps  = 0.5*k^2 ;    %imaginary shift of shifted Laplacian

%number of fine points and levels
%computed from number of points in 
%coarsest grid and wavenumber
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);  

k=0;

switch dim
    case 1
        A = helmholtz(k,0,npf,bc);
        M = helmholtz(k,eps,npf,bc);

    case 2
        A = helmholtz2(k,0,npf,npf,bc);
        M = helmholtz2(k,eps,npf,npf,bc);
end

% Test of V-cycle on a Helmholtz or Poisson problem

% Construct matrices on all grids and interpolation operators
op_type = 'gal';
[mg_mat,mg_split,restrict,interp] = mg_setup(A,k,0,op_type,numlev,bc,dim);
b  = zeros(length(A),1);
x0 = randn(length(A),1);


% Parameters of V-cycle and Jacobi iteration
npre = 1; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 10;
r0   = norm(b-A*x0);

%% Test of multigrid on Helmholtz problem
% Use the profiler to evaluate performance of code
profile on
[x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
profile off

r0 
r1 = norm(b-A*x_sol)/norm(b-A*x0) %Multigrid diverges on Helmholtz problems!

%pause;

%% Test of multigrid on shifted Laplacian problem
[mg_mat,mg_split,restrict,interp] = mg_setup(M,k,eps,op_type,numlev,bc,dim);

profile on
[x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
profile off

r0
r1 = norm(b-M*x_sol)/norm(b-M*x0)

