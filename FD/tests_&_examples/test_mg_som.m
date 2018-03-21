%%Test multigrid Sommerfeld BC's
clear all; 
close all;
clc;

k         = 100;                  %wavenumber
ppw       = 0.5;                  %min points per wavelength%
npcg      = 1;                    %number of points in coarsest grid
bc        = 'som1';               %boundary conditions
dim       = 2;                    %dimension
factoreps = 0.5;
poweps    = 2;
eps       = factoreps*k^poweps ;  %imaginary shift of shifted Laplacian

[npf,numlev] = fd_npc_to_npf(npcg,k,ppw);
op_type      = 'gal';

%Constructing the matrices and multigrid operators
profile on
[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npcg,numlev,bc,dim);

%right hand side and initial guess
A   = mg_mat{1}; ex_sol = ones(length(A),1); 
b   = A*ex_sol;
x0  = rand(length(A),1); 

%Parameters of multigrid solver
npre = 1; npos = 1; w = 0.7; smo = 'wjac'; numcycles = 30;

% %Running multigrid cycle for shifted Laplacian
 relres = zeros(numcycles+1,1); relres(1)=norm(b-A*x0); 
 relerr = zeros(numcycles+1,1); relerr(1)=norm(ex_sol);

 for i=1:numcycles
     [x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,1);
     relres(i+1)=norm(b-A*x_sol);
     x0=x_sol;
 end
 
 relres = relres/relres(1);
 factor = relres(2:length(relres))./relres(1:length(relres)-1);

 profile off

 relres
 
