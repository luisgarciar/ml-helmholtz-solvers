%% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian
% 2-D Example, Sommerfeld boundary conditions
clc
clear global;
close all;
npc = 3;    %number of interior points in coarsest grid in one dim
dim = 2; %boundary conditions, dimension

%wavenumber and imaginary shift of shifted Laplacian
k  = 40;
factoreps = 0.5; poweps = 2;
eps = factoreps*k^poweps; %Helmholtz problem

%% %Exact solution of Dirichlet problem
% u = @(x,y)  x.*(x-1).^2.*(y.^2).*(y-1);
% f = @(x,y) -k^2*x.*(y.^2).*((x-1).^2).*(y-1)...
%             -2*x.*((x-1).^2).*(3*y-1)...
%             -2*(y.^2).*(3*x-2).*(y-1);
%
% %Helmholtz and Shifted Laplacian Matrices and rhs
% dim = 2;
% M  = helmholtz2(k,eps,npf,npf,bc);
%
% npx = npf; npy=npf; hx = 1/(npx+1); hy = 1/(npy+1);
% xgrid = hx*(1:1:npx); ygrid = hy*(1:1:npy);
% np = npx*npy;
% [X,Y]= meshgrid(xgrid,ygrid);
%
% u_ex = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[np,1]); %exact solution
% b    = f(X,Y); b = b'; b = reshape(b,[np,1]);  %right hand side

%% Multigrid Setup
profile on
op_type = 'gal'; %type of coarse operators (galerkin or rediscretized)
[npf,numlev] = fem_npc_to_npf(npc,k,0.5);  %number of points in finest grid (1D)

%Construction of the linear system and the preconditioner
bc = 'som';
%Construct square mesh of meshsize h
h = 1/(npf+1);
[node,elem] = squaremesh([0,1,0,1],h);

%Find boundary nodes
[bdNode,bdEdge,isBdNode] = findboundary(elem);

%Sets Sommerfeld boundary conditions on all boundary edges
bdFlag = setboundary(node,elem,'ABC');

%The structures pde(helm,SL) contain data for the Helmholtz and
%shifted Laplace problems
pdehelm = helmholtz2Dconstantwndata(k,0,1);
pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);
option.tol = 1e-12;

%Helmholtz and shifted Laplace matrices

[mg_mat,mg_split,~,~]= mg_setupfem_2D(npc,numlev,pdeSL);
[mg_mat_helm,~,restrict,interp]= mg_setupfem_2D(npc,numlev,pdehelm);

A    = mg_mat_helm{1};
Aeps = mg_mat{1};
b    = ones(length(A),1);
u_gal = A\b; %exact solution

% Setting the CSL preconditioners Minv
% Parameters of MG-cycle and Jacobi iteration
b    = ones(length(A),1);
x0   = zeros(size(b));
npre = 1; npos = 1; w = 0.5; smo = 'wjac'; numcycles = 1;
Minv_wmg = @(v)feval(@Wcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,numcycles);
Minv_vmg = @(v)feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,2);
Minv_fmg = @(v)feval(@Fcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,numcycles);

% %Parameters of GMRES iteration
tol   = 1e-10;
maxit = 300;

% %GMRES iteration with LU preconditioner
[L,U] = lu(Aeps);
AMinv_lu = @(v) A*(U\(L\v));
tic
[xlu,flaglu,relreslu,iterlu,resveclu] = gmres(AMinv_lu,b,[],tol,maxit);
time_lu = toc;
u_lu =U\(L\xlu);

%GMRES iteration with V-cycle preconditioner
AMinv_vmg = @(v) A*feval(Minv_vmg,v);
tic
[xv,flagv,relresv,iterv,resvecv] = gmres(AMinv_vmg,b,[],tol,maxit);
time_v = toc;
u_vmg = Minv_vmg(xv);

% GMRES iteration with F-cycle preconditioner
AMinv_fmg = @(v) A*feval(Minv_fmg,v);
tic
[xf,flagf,relresf,iterf,resvecf] = gmres(AMinv_fmg,b,[],tol,maxit);
time_f = toc;
u_fmg = Minv_fmg(xf);

% GMRES iteration with W cycle preconditioner
AMinv_wmg =  @(v) A*feval(Minv_wmg,v);
tic
[xw,flagw,relresw,iterw,resvecw] = gmres(AMinv_wmg,b,[],tol,maxit);
time_w = toc;
u_wmg = Minv_wmg(xw);


%% Plots
figure(1)

semilogy(1:(iterlu(2)+1),resveclu'/resveclu(1),'r-+');
hold on
semilogy(1:(iterv(2)+1),resvecv'/resvecv(1),'r-+');
semilogy(1:(iterf(2)+1),resvecf'/resvecf(1),'k-*');
semilogy(1:(iterw(2)+1),resvecw'/resvecw(1),'b-*');

legend('LU preconditioner','Vcycle preconditioner','Fcycle preconditioner',...
    'Wcycle preconditioner');

time_v
time_f
time_w


norm(u_gal-u_wmg,Inf)
norm(u_gal-u_vmg,Inf)
norm(u_gal-u_fmg,Inf)




