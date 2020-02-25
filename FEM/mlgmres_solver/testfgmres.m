%% test fgmres

% n = 20;
% A1 = randn(n,1);
% A2 = randn(n,n);
% b  = ones(n,1);
% A  = A1+1i*A2;
% x1 = gmres(A,b);
% x2 = fgmres(A,b,1e-8);
% 
% res1 = norm(b-A*x1);
% res2 = norm(b-A*x2);

 
%% test with Helmholtz matrix and Shifted Laplace preconditioner

dim       = 2;
poweps    = 2;
factoreps = 1;

k  = 40;
bc = 'som';

reflevs   = 1;   %
restart   = [];
tol       = 1e-6;
maxit     = 200;

npcc = 3;
par  = 0.7;
npf  = ceil(k^(3/2));
np   = npf-2;

[npf,numlev] = fem_npc_to_npf(npcc,k,par);  %number of points in finest grid (1D)

h = 1/(npcc+1);
[node,elem] = squaremesh([0 1 0 1],h);

%refining the mesh numlev times
for j = 1:numlev-1
    [node,elem] = uniformrefine(node,elem);
end

%set boundary conditions
[bdNode,bdEdge,isBdNode] = findboundary(elem);
bdFlag = setboundary(node,elem,'ABC');

pdehelm    = helmholtz2Dconstantwndata(k,0,1);
pdeSL      = helmholtz2Dconstantwndata(k,factoreps,poweps);
option.tol = 1e-8;

[eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
[eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);

[mg_mat,mg_split,restr,interp]= mg_setupfem_2D(npcc,numlev,pdeSL);

A      = eqn1.A;
Aeps  = mg_mat{1};
%[L,U]  = lu(Aeps);

n1 = length(A);
n2 = length(Aeps);

assert(n1==n2,'incorrect matrix size')

b       = ones(length(A),1);
restart = [];
u0      = zeros(length(A),1);

npre = 1; npos = 1; w  = 0.7; numit = 20; smo = 'wjac';
Aepsinv = @(x) Fcycle(mg_mat,mg_split,restr,interp,u0,x,npre,npos,w,smo,1);
mat1 = @(x) A*Aepsinv(x);

tol = 1e-8;
Aa = @(x,tol) A*x;
P = @(x,tol) Aepsinv(x);

[x1,iter,resids] = fgmres(Aa, b, tol, 'P',P,'max_iters', 1,'restart',200);
[x2, ~, ~, iter1, ~] = gmres(mat1,  b, restart, tol, maxit);




