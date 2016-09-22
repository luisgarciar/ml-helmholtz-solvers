%Test multigrid Sommerfeld BC's
clc; clear all; close all;
bc  = 'dir';
dim = 2;
k   = 0;
eps = k^2*0.5;
npc = 5; 
ppw = 25;
[npf,lev] = fd_npc_to_npf_som(npc,k,ppw); 

%Constructing the matrices and multigrid operators
A = helmholtz2(k,0,npf,npf,bc);
%M = helmholtz2(k,eps,npf,npf,bc);

[grid_matrices,grid_smooth,restrict,interp] = mg_setup_som(A,lev,bc,dim);

%right hand side and initial guess
b   = randn(length(A),1);  b=b/norm(b);
ex_sol = A\b;

x0  = zeros(length(A),1);
%x0  = ex_sol; 

%Parameters of multigrid solver
npre = 1; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 20;

%Running multigrid cycle for shifted Laplacian
relres = zeros(numcycles+1,1); relres(1)=norm(b-A*x0);

for i=1:numcycles
    [x_sol] = Vcycle_som(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,1);
    relres(i+1)=norm(b-A*x_sol);
    x0=x_sol;
end

relres=relres/relres(1)

