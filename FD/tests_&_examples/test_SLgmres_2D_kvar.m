%% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian 
% 2-D Example, Dirichlet boundary conditions
clear all; close all;

npc = 3;    %number of interior points in coarsest grid in one dim 
bc  = 'dir'; dim = 2; %boundary conditions, dimension

%variable wavenumber and imaginary shift of shifted Laplacian
kref    = 10; 
kvar    = @(x,y) klay(x,y,kref);
epsvar  = @(x,y) 0.5*(klay(x,y,kref).^2); 
zero    = @(x,y) 0*x;
k       = 5*kref;
ppw     = 12;                          %number of points per wavelength
[npf,lev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)

%Grid for plotting
npx = npf; npy=npf; hx = 1/(npx+1); hy = 1/(npy+1);
xgrid = hx*(1:1:npx); ygrid = hy*(1:1:npy);
np    = npx*npy;
[X,Y] = meshgrid(xgrid,ygrid);

%Helmholtz and Shifted Laplacian Matrices
% Dirichlet BC
A        = helmholtz2var(kvar,zero,npx,npy,bc);
M        = helmholtz2var(kvar,epsvar,npx,npy,bc);
b        = zeros(length(A),1); 
center   = npx*floor(npy/2)+ceil(npx/2);
b(center)= 1/(hx*hy);
%u_ex = A\b; sol=reshape(real(u_ex),npx,npy);
%surf(X,Y,sol); title('Real part of solution')
%pause(3); close all
%return

%Helmholtz and Shifted Laplacian Matrices (only test)
% bc = 'som'
% xgrid = hx*(0:1:npx+1); ygrid = hy*(0:1:npy+1);
% np    = (npx+1)*(npy+1);
% [X,Y] = meshgrid(xgrid,ygrid);
% A     = helmholtz2var(kvar,zero,npx,npy,bc);
% M     = helmholtz2var(kvar,epsvar,npx,npy,bc);
% b     = zeros(length(A),1); 
% center = npx*floor(npy/2)+ceil(npx/2); b(center)= 1/(hx*hy);
% u_ex = A\b; sol=reshape(real(u_ex),(npx+2),(npy+2));
% surf(X,Y,sol); title('Real part of solution')
% return

%% Multigrid Setup
tic
[galerkin_matrices,galerkin_split,restrict,interp]= mg_setup(M,lev,bc,dim);
setuptime=toc;

%Setting the SL preconditioner Minv
%Parameters of V-cycle and Jacobi iteration
b    = zeros(length(A),1); b(ceil(length(A)/4),1)=1/(hx*hy);
x0   = zeros(length(A),1);
npre = 2; npos = 0; w = 2/3; smo = 'wjac'; numcycles = 1;
Minv = @(v)feval(@Vcycle,galerkin_matrices,galerkin_split,restrict,interp,x0,v,npre,npos,w,smo,numcycles);
 
%Parameters of GMRES iteration
tol   = 1e-10;
maxit = length(A);

%GMRES iteration without preconditioner (too slow for large wavenumbers)
% tic
% [~,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
% time1 = toc;
 
% GMRES iteration with left SL preconditioner
tic
[x2,flag2,relres2,iter2,resvec2] = gmres(A,b,[],tol,maxit,Minv);
time2 = toc;

% GMRES iteration with right SL preconditioner
AMinv = @(v) A*feval(Minv,v);
tic
[x3,flag3,relres3,iter3,resvec3] = gmres(AMinv,b,[],tol,maxit);
time3 = toc;

%semilogy(1:(iter1(2)+1),resvec1'/resvec1(1),'b-+');
%hold on
figure(1)
semilogy(1:(iter2(2)+1),resvec2'/resvec2(1),'r-+');
hold on
semilogy(1:(iter3(2)+1),resvec3'/resvec3(1),'k-*');

figure(2)
u_ex = A\b; 
sol = reshape(real(u_ex),npx,npy); surf(X,Y,sol); title('Real part of solution')
figure(3)
sol = reshape(real(x2),npx,npy);   surf(X,Y,sol); title('Real part of solution')
