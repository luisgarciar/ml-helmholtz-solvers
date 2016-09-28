%%Test multigrid Sommerfeld BC's
clc; clear all; close all;
bc  = 'dir';
dim = 2;
k   = 100;
eps = 0;
npc = 3; 
ppw = 20;

[npf,lev] = fd_npc_to_npf_som(npc,k,ppw); 

%Constructing the matrices and multigrid operators
k=0;
A = helmholtz2(k,eps,npf,npf,bc);
[galerkin_matrices,galerkin_splitting,restrict,interp] = mg_setup_som(A,lev,bc,dim);

figure(1);
spy(galerkin_matrices{1})

figure(2)
spy(galerkin_matrices{2})

%%
%right hand side and initial guess
b   = zeros(length(A),1); 
ex_sol = A\b;

x0  = rand(length(A),1); 
%x0  = ex_sol; 

%Parameters of multigrid solver
npre = 1; npos = 1; w = 0.6; smo = 'gs'; numcycles = 20;

%Running multigrid cycle for shifted Laplacian
relres = zeros(numcycles+1,1); relres(1)=norm(b-A*x0);

for i=1:numcycles
    [x_sol] = Vcycle_som(galerkin_matrices,galerkin_splitting,restrict,interp,x0,b,npre,npos,w,smo,1);
    relres(i+1)=norm(b-A*x_sol);
    x0=x_sol;
end

relres = relres/relres(1)
npf = 5;

factor = relres(2:length(relres))./relres(1:length(relres)-1)
figure(1);
spy(galerkin_matrices{1})

figure(2)
spy(galerkin_matrices{2})

P = perm_rb_som(length(A));