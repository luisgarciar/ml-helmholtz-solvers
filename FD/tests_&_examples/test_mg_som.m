%%Test multigrid Sommerfeld BC's
clear all; 
close all;
clc;

k    = 80;     %wavenumber
ppw  = 25;     %min points per wavelength%
npc  = 3;      %number of points in coarsest grid
bc   = 'dir';  %boundary conditions
dim  = 2;      %dimension
eps  = 0.5*k^2 ;   %imaginary shift of shifted Laplacian

[npf,numlev] = fd_npc_to_npf(npc,k,ppw); 
op_type = 'rd';

%k=0; eps=0;

%Constructing the matrices and multigrid operators
A = helmholtz2(k,eps,npf,npf,bc);
[mg_mat,mg_split,restrict,interp] = mg_setup(A,k,eps,op_type,numlev,bc,dim);

%figure(1);
%spy(mg_mat{1})

%figure(2);
%spy(mg_mat{2})

%%
%right hand side and initial guess
b      = ones(length(A),1); 
ex_sol = A\b;

x0  = rand(length(A),1); 
%x0  = ex_sol; 

%Parameters of multigrid solver
npre = 1; npos = 1; w = 0.6; smo = 'wjac'; numcycles = 8;

%Running multigrid cycle for shifted Laplacian
relres = zeros(numcycles+1,1); relres(1)=norm(b-A*x0); 
relerr = zeros(numcycles+1,1); relerr(1)=norm(ex_sol);

for i=1:numcycles
    [x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,1);
    relres(i+1)=norm(b-A*x_sol);
    x0=x_sol;
end

relres = relres/relres(1)
factor = relres(2:length(relres))./relres(1:length(relres)-1)


%%
%npf = 5;
%figure(1);
%spy(mg_mat{1})
%figure(2)
%spy(mg_mat{2})
%P = perm_rb(length(A));