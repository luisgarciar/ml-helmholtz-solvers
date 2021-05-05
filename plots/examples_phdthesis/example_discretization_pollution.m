dim       = 2;
poweps    = 2;
factoreps = 1;

npcc = 3;
par  = 0.7;
npf  = ceil(k^(3/2));
np   = npf-2;

[npf,numlev] = fem_npc_to_npf(npcc,k,par);  %number of points in finest grid (1D)

h = 1/(npcc+1);
[node,elem] = squaremesh([0 1 0 1],h);

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
  Aeps   = eqn.A;
 
%Right hand side
 b  = eqn1.A;
 x0 = zeros(size(b));
 
%Solving the Galerkin problem 
[u_gal,~,~,iter2] = mlfgmres(b,x0,mg_matHelm,mg_matCSL,mg_splitCSL,restr,interp,maxiter,tol);


%Computing the relative error of the Galerkin problem w.r.t.
%the exact solution

%We need the 






