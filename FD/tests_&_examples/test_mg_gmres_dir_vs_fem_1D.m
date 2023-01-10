%% Test multigrid and GMRES for 1D Problems,
% Comparison between Dirichlet (FD) and Mixed BCs (Som)
%
%Results:
%Using a constant (equal to one) function on the right hand side, 
%the Sommerfeld problem requires more iterations than the Dirichlet problem.
%to reduce the relative residual below a given tolerance.
%
%Using a point source on the right hand side, the Sommerfeld problem
%requires fewer iterations (as remarked by the experiments of several authors)

clear global; 
close all;
clc;

% Parameters
k        = 40;   %wavenumber
ppw      = 30;    %min points per wavelength%
npc_int  = 3;     %number of interior points in coarsest grid
dim      = 1;     %dimension

% Parameters of multigrid solver
npre = 2; npos = 1; w = 2/3; smo = 'wjac';

% Number of points in fine grid
[npf,numlev] = fd_npc_to_npf(npc_int,k,ppw);
% Number of points in fine grid for FEM problem
npc_fem = npc_int+1;
npf_fem = npf+1;

%% Mixed problem (Dirichlet + Sommerfeld), FEM
op_type      = 'gal'; %galerkin coarse operators

%Helmholtz FEM matrix
A_fem = helmholtzfem(k,npf_fem,0);

%Multigrid setup for shifted Laplacian
profile on
eps = 0.5*k^2;
[mg_mat_fem,mg_split_fem,restrict_fem,interp_fem] = mg_setupfem(k,eps,op_type,npc_fem,numlev,dim);

%right hand side and initial guess
M_fem  = mg_mat_fem{1};

%constant right hand side
%rhs_fem  = ones(npf_fem,1); rhs_fem(npf_fem) = 0.5; h = 1/npf_fem; rhs_fem = h*rhs_fem; %(constant function f=1)

%point source on right hand side
rhs_fem = ones(npf_fem,1); rhs_fem(npf_fem/2)=h*1;
x0_fem = zeros(size(rhs_fem));


%% Multigrid comparison

%Mixed Problem, FEM

%multigrid cycle for shifted Laplacian
numcycles = 12;
relres_fem_mg = zeros(numcycles+1,1); relres_fem_mg(1)=norm(rhs_fem-M_fem*x0_fem); 
%relerr_fem = zeros(numcycles+1,1); %relerr_fem(1)=norm(ex_sol_dir);

for i=1:numcycles
    [x_sol] = Vcyclefem(mg_mat_fem,mg_split_fem,restrict_fem,interp_fem,x0_fem,rhs_fem,npre,npos,w,smo,1);
    relres_fem_mg(i+1)=norm(rhs_fem-M_fem*x_sol);
    x0_fem =x_sol;
end
 
relres_fem_mg = relres_fem_mg/relres_fem_mg(1);
factor_fem = relres_fem_mg(2:length(relres_fem_mg))./relres_fem_mg(1:length(relres_fem_mg)-1);
profile off

% Dirichlet problem
bc = 'dir';

% Multigrid setup for shifted Laplacian
[mg_mat_dir,mg_split_dir,restrict_dir,interp_dir] = mg_setup(k,eps,op_type,npc_int,numlev,bc,dim);

% Matrix, right hand side and initial guess
M_dir =  mg_mat_dir{1};

% constant right hand side
%rhs_dir =  ones(length(M_dir),1); 

% point source
rhs_dir =  zeros(length(M_dir),1);  rhs_dir(floor(npf/2))=1;
x0 =  rand(length(M_dir),1); 

% Setup of residuals and error
relres_mg_dir = zeros(numcycles+1,1); relres_mg_dir(1) = norm(rhs_dir-M_dir*x0); 
%relerr_dir = zeros(numcycles+1,1); relerr_dir(1) = norm(ex_sol_dir);

% Multigrid cycle of Dirichlet problem
tic
for i=1:numcycles
    [x_sol] = Vcycle(mg_mat_dir,mg_split_dir,restrict_dir,interp_dir,x0,rhs_dir,npre,npos,w,smo,1);
    relres_mg_dir(i+1)=norm(rhs_dir-M_dir*x_sol);
    x0 =x_sol;
end
time_dir = toc;  

relres_mg_dir = relres_mg_dir/relres_mg_dir(1);
factor_dir = relres_mg_dir(2:length(relres_mg_dir))./relres_mg_dir(1:length(relres_mg_dir)-1);

relres_fem_mg
relres_mg_dir

factor_fem
factor_dir


%% Test of MG + GMRES (Shifted Laplacian + Helmholtz)
% Dirichlet problem
dim = 1;
A_dir  = helmholtz(k,0,npf,bc);  %set eps=0 for Helmholtz problem
%M  = helmholtz(k,eps,npf,bc); %shifted Laplacian
 
%% Multigrid Setup
% Construct matrices on all grids and interpolation operators
op_type = 'gal';
[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc_int,numlev,bc,dim);

M = mg_mat{1};
x0 = zeros(length(A_dir),1);
 

%% Test of Preconditioned GMRES
% Setting the MG preconditioner Minv
% Parameters of V-cycle and Jacobi iteration
npre = 2; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 1;
Minv_dir = @(v) feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,1);
AMinv_dir = @(v) A_dir*feval(Minv_dir,v);

% Parameters of GMRES iteration
tol   = 1e-10;
maxit = min(150,length(A_dir));
 
% GMRES iteration with left SL preconditioner
tic
[x_lgdir,flag_lgdir,relres_lgdir,iter_lgdir,resvec_lgdir] = gmres(A_dir,rhs_dir,[],tol,maxit,Minv_dir);
time_lgdir = toc;

% GMRES iteration with right SL preconditioner
tic
[x_rgdir,flag_rgdir,relres_rgdir,iter_rgdir,resvec_rgdir] = gmres(AMinv_dir,rhs_dir,[],tol,maxit);
time_rgdir = toc;

% Mixed BC Problem
x0 = zeros(length(A_fem),1); 
numcycles = 1;
Minv_fem = @(v) feval(@Vcycle,mg_mat_fem,mg_split_fem,restrict_fem,interp_fem,x0,v,npre,npos,w,smo,numcycles);
AMinv_fem = @(v) A_fem*feval(Minv_fem,v);
 
% GMRes parameters
tol   = 1e-10;
maxit = min(150,size(A_fem,1));

% GMRES + left preconditioning with MG - Shifted Laplacian  
tic
[x_lgfem,flag_lgfem,relres_lgfem,iter_lgfem,resvec_lgfem] = gmres(A_fem,rhs_fem,[],tol,maxit,Minv_fem);
time_lgfem = toc;

tic
[x_rgfem,flag_rgfem,relres_rgfem,iter_rgfem,resvec_rgfem] = gmres(AMinv_fem,rhs_fem,[],tol,maxit);
time_rgfem = toc;
 
iter_lgfem
iter_lgdir

iter_rgfem
iter_rgdir

time_lgfem
time_rgfem

time_lgdir
time_rgdir


%% Plots of GMRES results
figure(1)
semilogy(0:iter_rgfem(2), resvec_rgfem/resvec_rgfem(1), 'k-');
hold on
semilogy(0:(length(resvec_rgdir)-1), resvec_rgdir/resvec_rgdir(1), 'b-');
ylabel('relative residual')
xlabel('iteration')
 
legend('FEM problem with mixed BCs', 'FD problem with Dirichlet BCs')
title(['Solving the Helmholtz 1D equation with GMRES and CSL-preconditioner (k=',num2str(k),')'])
