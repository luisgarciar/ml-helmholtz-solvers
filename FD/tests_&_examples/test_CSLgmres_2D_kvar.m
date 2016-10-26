%% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian 
% 2-D Example, Dirichlet boundary conditions
clc
clear global;
close all;

npc = 1;    %number of interior points in coarsest grid in one dim 
bc  = 'dir'; dim = 2; %boundary conditions, dimension

%variable wavenumber and imaginary shift of shifted Laplacian
kref = 40; eps = 0.5;
ppw  = 20;  %number of points per wavelength
[npf,numlev] = fd_npc_to_npf(npc,kref,ppw);  %number of points in finest grid (1D)

kvar    = @(x,y) klay(x,y,kref);
epsvar  = @(x,y) eps*(klay(x,y,kref).^2); 
zero    = @(x,y) 0*x;

%Grid for plotting
npx   = npf; npy=npf; hx=1/(npx+1); hy=1/(npy+1);
xgrid = hx*(1:1:npx); ygrid = hy*(1:1:npy);
np    = npx*npy;
[X,Y] = meshgrid(xgrid,ygrid);

%Helmholtz and Shifted Laplacian Matrices
% Dirichlet BC
A   = helmholtz2var(kvar,zero,npx,npy,bc);
M   = helmholtz2var(kvar,epsvar,npx,npy,bc);
%b(index)= 1/(hx*hy);
%u_ex = A\b; sol = reshape(real(u_ex),npx,npy);
%surf(X,Y,sol); title('Real part of solution');

%% Multigrid Setup
tic
op_type   = 'gal';
[mg_mat,mg_split,restrict,interp] = mg_setup_kvar(kref,eps,op_type,npc,numlev,bc,dim);
setuptime = toc;

%Setting the SL preconditioner Minv
%Parameters of V-cycle and Jacobi iteration
b    = zeros(length(A),1); b(ceil(length(A)/2),1) = 1/(hx*hy);
x0   = zeros(length(A),1);
npre = 1; npos = 1; w = 1/3; smo ='wjac'; numcycles=2;

[L,U] = lu(M);
Minv_mg = @(v)feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,numcycles);
Minv_lu = @(v) L\(U\v);

%Parameters of GMRES iteration
tol   = 1e-6;
maxit = length(A);

%GMRES iteration without preconditioner (too slow for large wavenumbers)
% tic
% [~,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
% time1 = toc;
 
% GMRES iteration with left SL preconditioner
tic
[x2,flag2,relres2,iter2,resvec2] = gmres(A,b,[],tol,maxit,Minv_mg);
time2 = toc;

% GMRES iteration with right SL preconditioner
AMinv_mg = @(v) A*feval(Minv_mg,v);
tic
[x_rmg,flag_rmg,relres_rmg,iter_rmg,resvec_rmg] = gmres(AMinv_mg,b,[],tol,maxit);
time_rmg = toc;

AMinv_lu = @(v) A*feval(Minv_lu,v);
tic
[x_rlu,flag_rlu,relres_rlu,iter_rlu,resvec_rlu] = gmres(AMinv_lu,b,[],tol,maxit);
time_rlu = toc;

%semilogy(1:(iter1(2)+1),resvec1'/resvec1(1),'b-+');
%hold on
figure(1)
semilogy(1:(iter_rmg(2)+1),resvec_rmg'/resvec_rmg(1),'r-+');
hold on
semilogy(1:(iter_rlu(2)+1),resvec_rlu'/resvec_rlu(1),'k-*');



%figure(2)
%u_ex = A\b; 
%sol = reshape(real(u_ex),npx,npy); surf(X,Y,sol); title('Real part of solution')
figure(3)
sol = reshape(real(x2),npx,npy);   surf(X,Y,sol); title('Real part of solution')
