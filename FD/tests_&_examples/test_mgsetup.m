%Test of functions mg_setup, lin_interpol, fwrestriction
clear all; 
close all;
clc;

k    = 0;           %wavenumber
ppw  = 10;          %min points per wavelength%
npc  = 3;            %number of points in coarsest grid
bc   = 'dir';       %boundary conditions
dim  =  2;          %dimension
eps  = 0.5*k^2 ;  %imaginary shift of shifted Laplacian

%number of fine points and levels
%computed from number of points in 
%coarsest grid and wavenumber
[npf,lev] = fd_npc_to_npf(npc,k,ppw);  

switch dim
    case 1
        A = helmholtz(k,npf,bc);
        M = shift_laplace(k,b1,b2,npf,bc);

    case 2
        A = helmholtz2(k,npf,npf,bc);
        M = shift_laplace2(k,b1,b2,npf,npf,bc);
end

% Test of V-cycle on a Helmholtz or Poisson problem

% Construct matrices on all grids and interpolation operators
[grid_matrices,grid_smooth,restrict,interp] = mg_setup(A,lev,bc,dim);
 b  = zeros(length(A),1);
 x0 = randn(length(A),1);

% Parameters of V-cycle and Jacobi iteration
npre = 3; npos = 3; w = 2/3; smo = 'wjac'; numcycles = 10;
r0   = norm(b-A*x0);


%% Test of multigrid on Helmholtz problem
% Use the profiler to evaluate performance of code
profile on
[x_sol] = Vcycle(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
profile off

r0 
r1 = norm(b-A*x_sol)/norm(b-A*x0) %Multigrid diverges on Helmholtz problems!

pause;

%% Test of multigrid on shifted Laplacian problem
[grid_matrices,grid_smooth,restrict,interp] = mg_setup(M,lev,bc,dim);

profile on
[x_sol] = Vcycle(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
profile off

r0
r1 = norm(b-A*x_sol)/norm(b-A*x0);

