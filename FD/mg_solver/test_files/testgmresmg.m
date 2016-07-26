%% Test: Solving the Helmholtz eqn. with  GMRES preconditioned by multigrid and the shifted Laplacian

%% Wavenumber and gridsize data from the paper 
%Erlangga, Osterlee, Vuik, Appl. Num. Math., 2006
% kk   = [10 20 30 40 50]'
% npp = [31 63 95 127 192]' 
%lev = ceil(log2(k*ppw/2*pi)); npf = 2^lev-1; %number of 1D interior gridpoints
%lev = lev;

%% Run this section for 1D test
clear all; close all; clc
bc   = 'dir'; %boundary conditions ('som' to be implemented)
k    =  50;
ppw  =  20;  %min number of points per wavelength
b1   =  1;   b2 = 0.5;
npcc = 5;   %number of points on coarsest grid
[npf,lev] = npcc2npf(npcc,k,ppw);

%Helmholtz and Shifted Laplacian Matrices and right hand side
dim = 1;
A = helmholtz(k,npf,bc);
M = shift_laplace(k,b1,b2,npf,bc);
b = zeros(length(A),1); b(ceil(length(A)/2),1)=1;
z = zeros(length(A),1);

%Creating the matrices and intergrid operators on all levels
[SLgrid_matrices,SLgrid_smooth,restrict,interp] = mgsm_setup(M,lev,bc,dim);

%Multigrid Parameters
npre = 3; npos = 3; w  = 2/3; smo = 'gs'; numcycles=1;

%Setting the Shifted Laplace preconditioner Minv
Minv = @(x)feval(@Vcycle,SLgrid_matrices,SLgrid_smooth,...
                 restrict,interp,z,x,npre,npos,w,smo,numcycles); 

%Setting the Preconditioned Matrix AM^{-1}
AMinv = @(x)A*feval(Minv,x); 
                       
% Parameters of GMRES iteration
tol   = 1e-8;
maxit = length(A);
 
% GMRES iteration without preconditioner 
tic
[~,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
time_noprec=toc;

%GMRES iteration with SL preconditioner
tic 
[x,flag2,relres2,iter2,resvec2] = gmres(AMinv,b,[],tol,maxit);
time_prec=toc;

%Plots of the residuals
semilogy(1:(max(iter1)+1),resvec1/resvec1(1),'k--',...
         1:(max(iter2)+1),resvec2/resvec2(1),'r-+');
     
time_noprec
time_prec
     
%% 2D test

% Helmholtz and Shifted Laplacian Matrices and rhs
dim = 2;
k   = kk(4);
npf = npp(4);

A = helmholtz2(k,npf,npf,bc);
M = shift_laplace2(k,b1,b2,npf,npf,bc);
b = zeros(length(A),1); b(ceil(length(A)/2),1)=1;
z = zeros(length(A),1);

% Creating the matrices and intergrid operators on all levels
[SLgrid_matrices,SLgrid_smooth,restrict,interp] = mgsm_setup(M,lev,bc,dim);

% Multigrid Parameters
npre = 2; npos = 2; w  = 2/3; smo = 'wjac'; numcycles=1;

%Setting the SL preconditioner Minv
Minv = @(x)feval(@Vcycle2,SLgrid_matrices,SLgrid_smooth,...
                restrict,interp,z,x,npre,npos,w,smo,numcycles); 

%Setting the Preconditioned Matrix
AMinv = @(x)A*feval(Minv,x) 
                       
% Parameters of GMRES iteration
 tol   = 1e-6;
 maxit = length(A);
 
% GMRES iteration without preconditioner

%  tic
%  [~,flag1,relres1,iter1,resvec1] = ...
%    gmres(A,b,[],tol,maxit,[]);
%  time_noprec2d=toc;

%GMRES iteration with SL preconditioner
tic 
[x,flag2,relres2,iter2,resvec2] = ...
   gmres(AMinv,b,[],tol,maxit);
time_prec2d=toc;

%Plots of the residuals
% semilogy(1:(max(iter1)+1),resvec1/resvec1(1),'b-+',...
%          1:(max(iter2)+1),resvec2/resvec2(1),'r-+');
   
 semilogy(1:(max(iter2)+1),resvec2/resvec2(1),'r-+');
 
%time_noprec2d
time_prec2d

     


