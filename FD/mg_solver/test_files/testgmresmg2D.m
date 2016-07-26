%% Test of GMRES preconditioned by multigrid and the shifted Laplacian

%% Test in 2D
npcc = 5;                     %number of points in coarsest grid (1D) 
k    = 20; 
[npf,lev]  = npcc2npf(npcc,k,ppw);  %number of points in finest grid (1D)
bc   = 'dir';
b1   = 1; b2 = 0.5;
k   =  0;

%% Helmholtz and Shifted Laplacian Matrices and rhs
dim = 2;
A   = helmholtz2(k,npf,npf,bc);
M   = shift_laplace2(k,b1,b2,npf,npf,bc);
b   = zeros(length(A),1); b(ceil(length(A)/4),1)=1;
z   = zeros(length(A),1); 

%% Multigrid Setup
 tic
 [SLgrid_matrices,restrict,interp] = mgsm_setup(M,lev,bc,dim);
 setuptime = toc;

%%
%Multigrid Parameters
nu = 2; mu = 2; w  = 0.4;  %2/3

%Setting the SL preconditioner Minv
%Minv = @(v)feval(@Vcycle,SLgrid_matrices,restrict,interp,z,v,nu,mu,w);
[L,U] = lu(M);
Minv  = @(v) U\(L\v);

%Parameters of GMRES iteration
tol   = 1e-8;
maxit = length(A);
 
%GMRES iteration without preconditioner
% tic
% [~,flag1,relres1,iter1,resvec1] = ...
%   gmres(A,b,[],tol,maxit,[]);
% time1 = toc;

%GMRES iteration with M\x (LU) as preconditioner


%GMRES iteration with SL preconditioner
tic
[x,flag2,relres2,iter2,resvec2] = ...
   gmres(A,b,[],tol,maxit,Minv);
time3 = toc;
 
%semilogy(1:(max(iter2)+1),resvec2/resvec2(1),'r-+',1:(max(iter1)+1),resvec1/resvec1(1),'b-+');
