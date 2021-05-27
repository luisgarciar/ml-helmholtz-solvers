clear all;

dim       = 2;
poweps    = 2;
factoreps = 10;
k = 20;

npcc = 3;
par  = 0.7;
npf  = ceil(k^(3/2));
np   = npf-2;

[npf,numlev] = fem_npc_to_npf(npcc,k,par);  %number of points in finest grid (1D)

h = 1/(npf+1);
[node,elem] = squaremesh([0 1 0 1],h);
 
% %Refining the grid
% for j=1:numlev-1
% [node,elem] = uniformrefine(node,elem);
% end

[bdNode,bdEdge,isBdNode] = findboundary(elem);
bdFlag = setboundary(node,elem,'ABC');

% This section creates the data for the Helmholtz problem with exact solution  
%  -div(grad u)-k^2 u = 0 in Omega= (0,1)x(0,1)
%  grad(u) dot n - i*ku = g in bd(Omega) 
%  with a plane wave as exact solution, i.e.,
%  u = e^{i kk \cdot (x,y)}  where kk = k (cos(t), sin(t))
%  and t is the angle (direction) of the wave
%   Output: 
%   pde:    struct containing the following data:
%   All function handles to be applied to input of size (N,2)
%   'f':    function handle for right hand side (equal to 0)
%   'exactu': function handle for exact solution
%   'gradu': function handle for gradient of exact solution
%   'k':  wavenumber
%   'g': function handle for boundary data

t = pi/2; % direction of propagation
pdehelm = helmholtz2Dplanewavedata(k,t);
u_ex    = pdehelm.exactu;
f = pdehelm.f;
g = pdehelm.g; 

% Data for shifted Laplace problem (no need for right hand side)
pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);
option.tol = 1e-8;

% Assembly of the matrix and rhs vector
% stored in structs eqn1, eqn2
[eqn1,~] = helmholtz2Dfem(node,elem,pdehelm,bdFlag,bdEdge);
[eqn2,~] = helmholtz2Dfem(node,elem,pdeSL,bdFlag,bdEdge);

%Setup of multilevel Krylov solver and shifted Laplacian
 [mg_matHelm,mg_splitHelm,restr,interp]    = mg_setupfem_2D(npcc,numlev,pdehelm);
 [mg_matCSL,mg_splitCSL,restrCSL,interpCSL]= mg_setupfem_2D(npcc,numlev,pdeSL);
 
 %Helmholtz Matrix and shifted Laplace matrix
  A      = eqn1.A;
  Aeps   = eqn2.A;
 
%Right hand side
 b  = eqn1.b;
 
 %b = sparse(length(b),1); 
 %b(ceil(length(b)/2))=1;
 %x0 = sparse(length(b),1);
 
 maxiter = ones(length(mg_matHelm),1);
 maxiter(1:5)=[200,8,4,2,2]';
 tol = 1e-12;
 u_gal = A\b;
 
 
%Solving the Galerkin problem with the shifted Laplacian
npre = 0; npos = 1; w  = 0.6; smo = 'gs';
u0      = zeros(length(A),1);
Aepsinv = @(x,tol,num) Vcycle(mg_matCSL,mg_splitCSL,restrCSL,interpCSL,u0,x,npre,npos,w,smo,1);
AP  = @(x) A*Aepsinv(x,tol,10);

[L,U] = lu(Aeps);
Aepsinv2 = @(x,tol) U\(L\x);
AP2 =  @(x) A*Aepsinv2(x,tol);


tol = 1e-12;
maxit = 150;
restart = [];
print_level = 1;


%[u_iter0,flag1,relres1, iter1] = gmres(AP2,b,restart,tol,maxit);
%u_iter1 = Aepsinv2(u_iter0,tol);

%[u_iter2,flag2,relres2, iter2] = gmres(AP,b,restart,tol,maxit);
%u_iter3 = Aepsinv(u_iter2,tol,10);

[u_iter4,flag2,relres2, iter2] = gmres(A,b,restart,tol,maxit,Aepsinv);



%[u_iter1, iter, resid] = fgmres(A, b, tol, 'restart',200,'tol_exit', tol,'P',Aepsinv);

%Solving the Galerkin problem with mlpfgmres
%[u_iter2,~,~,iter2] = mlfgmres(b,u0,mg_matHelm,mg_matCSL,mg_splitCSL,restr,interp,maxiter,tol);

%plotting the real part of the solution
%showsolution(node,elem,real(full(u_gal)));
%showsolution(node,elem,real(full(u_iter2)));

% figure(1)
% showsolution(node,elem,real(full(u_iter1)));
% 
% figure(2)
% showsolution(node,elem,real(full(u_iter3)));


figure(2)
showsolution(node,elem,real(full(u_iter4)));

%Computing the relative error of the Galerkin problem w.r.t.
%the exact solution

%We need the 






