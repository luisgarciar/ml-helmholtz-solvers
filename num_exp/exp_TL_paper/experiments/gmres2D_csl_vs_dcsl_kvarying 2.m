%Comparison of Shifted Laplacian and two-level deflation for varying k
%Homogeneous problem

dim       = 2;
poweps    = 2;
factoreps = 1;
kk = [10 20 40];
m  = length(kk);

itercsl = zeros(m,1);
iterdef = zeros(m,1);


reflevs   = 1;   
restart   = [];
tol       = 1e-6;
maxit     = 200;
npcc = 3;
par  = 0.7;


for i = 1:m
    
k    = kk(i);
bc   = 'som';
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
Aeps   = mg_mat{1};

n1 = length(A);
n2 = length(Aeps);

assert(n1==n2,'incorrect matrix size')

b       = ones(length(A),1);
restart = [];
u0      = zeros(length(A),1);

npre = 1; npos = 1; w  = 0.7; numit = 20; smo = 'wjac';
Aepsinv = @(x) Fcycle(mg_mat,mg_split,restr,interp,u0,x,npre,npos,w,smo,1);
mat1 = @(x) A*Aepsinv(x);

[~, ~, ~, itercsl(i), ~] = gmres(mat1,  b, restart, tol, maxit);

R = restr{1};
P = interp{1};

Ac = R*A*P;

[Lc,Uc] = lu(Ac);
  
cgc =  @(x) P*(Uc\(Lc\(R*x)));
def =  @(x) Aepsinv(x-A*cgc(x))+cgc(x);
mat2 = @(x) A*def(x);
 
[~, ~, ~, iterdef(i), ~] = gmres(mat2,  b, restart, tol, maxit);

end