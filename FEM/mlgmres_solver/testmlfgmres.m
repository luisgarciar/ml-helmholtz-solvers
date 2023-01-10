%% test mlfgmres.m 
%% test with Helmholtz matrix and Shifted Laplace preconditioner
dim       = 2;
poweps    = 2;
factoreps = 1;

k  = 160;
bc = 'som';

maxit = 100;
npcc  = 1;
par   = 0.6;
npf   = ceil(k^(3/2));
np    = npf-2;

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

%building the matrices
[eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
[eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);

%setting up the multilevel structure
[mg_matHelm,mg_splitHelm,restr,interp]= mg_setupfem_2D(npcc,numlev,pdehelm);
[mg_matCSL,mg_splitCSL,restrCSL,interpCSL]= mg_setupfem_2D(npcc,numlev,pdeSL);

A      = mg_matHelm{1};
Aeps   = mg_matCSL{1};

n1 = length(A);
n2 = length(Aeps);

assert(n1==n2,'incorrect matrix size')

b       = ones(length(A),1);
restart = [];
u0      = zeros(length(A),1);

npre = 1; npos = 1; w  = 2/3; smo = 'wjac';
Aepsinv = @(x) Vcycle(mg_matCSL,mg_splitCSL,restrCSL,interpCSL,u0,x,npre,npos,w,smo,1);
mat1 = @(x) A*Aepsinv(x);

tol = 1e-6;
Aa  = @(x,tol) A*x;
AP  = @(x) A*Aepsinv(x);

maxiter = ones(length(mg_matHelm),1);
maxiter(1:5)=[20,6,2,2,1]';
x0 = zeros(n1,1);
 
tol = 1e-6;

  
tic
[x1,flag1,relres,iter]  = gmres(AP, b, [], tol, maxit);
timecsl = toc;

tic
[x2,flag2,resvec,iter2] = mlfgmres(b,x0,mg_matHelm,mg_matCSL,mg_splitCSL,restr,interp,maxiter,tol);
timeml = toc;
