%% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian 
% 2-D Example, Dirichlet boundary conditions
clc
clear global;
close all;

npc = 1;    %number of interior points in coarsest grid in one dim 
bc  = 'som'; dim = 2; %boundary conditions, dimension

%variable wavenumber and imaginary shift of shifted Laplacian
kref    = 80; eps = 0.5;
kvar    = @(x,y) klay(x,y,kref);
epsvar  = @(x,y) eps*(klay(x,y,kref).^2); 
zero    = @(x,y) 0*x;
ppw     = 12;                          %number of points per wavelength
[npf,lev] = fd_npc_to_npf(npc,kref,ppw);  %number of points in finest grid (1D)

%Grid for plotting
npx   = npf; npy=npf; hx=1/(npx+1); hy=1/(npy+1);
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
[galerkin_matrices,galerkin_split,restrict,interp]... 
= mg_setup_kvar(kref,eps,'gal',npc,lev,bc,dim);
setup_time=toc;

%Setting the SL preconditioner Minv
%Parameters of V-cycle and Jacobi iteration
b    = zeros(length(A),1); b(ceil(length(A)/2),1) = 1/(hx*hy);
x0   = zeros(length(A),1);
npre = 1; npos = 1; w = 0.6; smo = 'wjac'; numcycles = 1;
Minv_fmg = @(v)feval(@Fcycle,galerkin_matrices,galerkin_split,restrict,...
                            interp,x0,v,npre,npos,w,smo,numcycles);

Minv_vmg = @(v)feval(@Vcycle,galerkin_matrices,galerkin_split,restrict,...
                            interp,x0,v,npre,npos,w,smo,numcycles);

Minv_wmg = @(v)feval(@Wcycle,galerkin_matrices,galerkin_split,restrict,...
                            interp,x0,v,npre,npos,w,smo,numcycles);
                                                
%Parameters of GMRES iteration
tol   = 1e-8;
maxit = 300;

%GMRES iteration without preconditioner (too slow for large wavenumbers)
% tic
% [~,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
% time1 = toc;
 
% GMRES iteration with V cycle preconditioner
AMinv_vmg = @(v) A*feval(Minv_vmg,v);
tic
[xv,flagv,relresv,iterv,resvecv] = gmres(AMinv_vmg,b,[],tol,maxit);
time_v = toc;

% GMRES iteration with F cycle preconditioner
AMinv_fmg = @(v) A*feval(Minv_fmg,v);
tic
[xf,flagf,relresf,iterf,resvecf] = gmres(AMinv_fmg,b,[],tol,maxit);
time_f = toc;

% GMRES iteration with W cycle preconditioner
AMinv_wmg =  @(v) A*feval(Minv_wmg,v);
tic
[xw,flagw,relresw,iterw,resvecw] = gmres(AMinv_wmg,b,[],tol,maxit);
time_w = toc;

%semilogy(1:(iter1(2)+1),resvec1'/resvec1(1),'b-+');
%hold on
figure(1)
semilogy(1:(iterv(2)+1),resvecv'/resvecv(1),'r-+');
hold on
semilogy(1:(iterf(2)+1),resvecf'/resvecf(1),'k-*');
semilogy(1:(iterw(2)+1),resvecw'/resvecw(1),'b-*');

legend('Vcycle preconditioner','Fcycle preconditioner',...
        'Wcycle preconditioner');

time_v
time_f
time_w


% figure(2)
% u_ex = A\b; 
% sol = reshape(real(u_ex),npx,npy); surf(X,Y,sol); title('Real part of solution')
% figure(3)
% sol = reshape(real(x2),npx,npy);   surf(X,Y,sol); title('Real part of solution')
