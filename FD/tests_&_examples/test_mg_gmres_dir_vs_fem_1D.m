%% Test multigrid and GMRES, 
% 1D Problems, Dirichlet (FD) and Mixed BCs (Som)
clear global; 
close all;
clc;

%% Parameters
k     = 10;       %wavenumber
ppw   = 30;       %min points per wavelength%
npcg  = 5;        %number of points in coarsest grid
dim   = 1;        %dimension
eps   = 0.5*k^2 ; %imaginary shift of shifted Laplacian

% Parameters of multigrid solver
npre = 2; npos = 1; w = 2/3; smo = 'wjac';

%% Mixed problem (Dirichlet + Sommerfeld), FEM
bc           = 'som';      %boundary conditions
[npf,numlev] = fem_npc_to_npf(npcg,k,ppw);
op_type      = 'gal';

%Helmholtz FEM matrix
eps = 0;
A_fem = helmholtzfem(k,npf,eps);

%Multigrid setup for shifted Laplacian
profile on
eps = 0.5*k^2;
[mg_mat_fem,mg_split_fem,restrict_fem,interp_fem] = mg_setupfem(k,eps,op_type,npcg,numlev,dim);

%%
% right hand side and initial guess
M_fem = mg_mat_fem{1}; 
b_fem = ones(npf,1); b_fem(npf) = 0.5; h = 1/npf; b_fem = h*b_fem; %(constant function f=1)
x0_fem = zeros(size(b_fem));

numcycles = 12;
% multigrid cycle for shifted Laplacian Sommerfeld
relres_fem = zeros(numcycles+1,1); relres_fem(1)=norm(b_fem-M_fem*x0_fem); 
relerr_fem = zeros(numcycles+1,1); %relerr_fem(1)=norm(ex_sol_dir);

for i=1:numcycles
    [x_sol] = Vcycle(mg_mat_fem,mg_split_fem,restrict_fem,interp_fem,x0_fem,b_fem,npre,npos,w,smo,1);
    relres_fem(i+1)=norm(b_fem-M_fem*x_sol);
    x0_fem =x_sol;
end
 
relres_fem = relres_fem/relres_fem(1);
factor_fem = relres_fem(2:length(relres_fem))./relres_fem(1:length(relres_fem)-1);
profile off


%% Dirichlet problem
bc = 'dir';
eps = 0;
npcg_dir = npcg-1;
npf_dir = npc_numlev_to_npf(npcg_dir,numlev);
A_dir = helmholtz(k,eps,npf_dir,bc);

% Multigrid setup for shifted Laplacian
eps = 0.5*k^2;
[mg_mat_dir,mg_split_dir,restrict_dir,interp_dir] = mg_setup(k,eps,op_type,npcg_dir,numlev,bc,dim);

%Matrix, right hand side and initial guess
M_dir =  mg_mat_dir{1}; 
b_dir =  ones(length(M_dir),1);
ex_sol_dir =  A_dir\b_dir;
x0 =  rand(length(M_dir),1); 

%Setup of residuals and error
relres_dir = zeros(numcycles+1,1); relres_dir(1) = norm(b_dir-M_dir*x0); 
relerr_dir = zeros(numcycles+1,1); relerr_dir(1) = norm(ex_sol_dir);

tic
for i=1:numcycles
    [x_sol] = Vcycle(mg_mat_dir,mg_split_dir,restrict_dir,interp_dir,x0,b_dir,npre,npos,w,smo,1);
    relres_dir(i+1)=norm(b_dir-M_dir*x_sol);
    x0 =x_sol;
end
time_dir = toc;  

relres_dir = relres_dir/relres_dir(1);
factor_dir = relres_dir(2:length(relres_dir))./relres_dir(1:length(relres_dir)-1);

relres_fem
relres_dir

factor_fem
factor_dir


%% Test of MG + GMRES (Shifted Laplacian + Helmholtz)

%% Dirichlet problem
numcycles = 1;
Minv_dir   = @(v)feval(@Vcycle,mg_mat_dir,mg_split_dir,restrict_dir,interp_dir,x0,v,npre,npos,w,smo,numcycles);

%GMRes parameters
tol   = 1e-5;
maxit = min(150,length(A_dir));

tic
[x_gdir,flag_gdir,relres_gdir,iter_gdir,resvec_gdir] = gmres(A_dir,b_dir,[],tol,maxit,Minv_dir);
time_gdir = toc;


%% Mixed BC Problem
bc = 'som';
x0 = zeros(length(A_fem),1); 
numcycles = 3;
Minv_fem = @(v) feval(@Vcycle,mg_mat_fem,mg_split_fem,restrict_fem,interp_fem,x0,v,npre,npos,w,smo,numcycles);

%GMRes parameters
tol   = 1e-7;
maxit = min(150,size(A_fem,1));

tic
[x_gfem,flag_gfem,relres_gfem,iter_gfem,resvec_gfem] = gmres(A_fem,b_fem,[],tol,maxit,Minv_fem);
time_gfem = toc;

iter_gfem
iter_gdir


%% test of BiCGSTAB
%[~,FLAG,RELRES,ITER,RESVEC] = bicgstab(A_som,b_som,tol,maxit,Minv_fem);
%ITER


%% Plots of GMRES results
figure(1)
semilogy(0:iter_gfem(2), resvec_gfem/resvec_gfem(1), 'k-');
hold on
semilogy(0:(length(resvec_gdir)-1), resvec_gdir/resvec_gdir(1), 'b-');
ylabel('relative residual')
xlabel('iteration')

legend('Sommerfeld BCs', 'Dirichlet BCs')
title(['1D Helmholtz with CSL-preconditioner (k=',num2str(k),')'])

