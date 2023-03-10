% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian 
% 2-D Example, Sommerfeld boundary conditions
clc
clear global;
close all;

npc = 3;     %number of interior points in coarsest grid in one dim 
bc  ='som';  %boundary conditions
dim = 2;     %dimension

%wavenumber and imaginary shift of shifted Laplacian
k   = 40;  eps = 0.5*k^2; %Helmholtz problem
ppw = 20;                 %number of points per wavelength
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)


%% Multigrid Setup

A   = helmholtz2(k,0,npf,npf,bc); %Helmholtz matrix
tic
op_type = 'gal'; %galerkin operators
%Setting up the MG Operators
[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
setuptime = toc;

%Setting the SL preconditioner Minv
%Parameters of V-cycle and Jacobi iteration
b  = ones(length(A),1);
x0 = zeros(size(b));
npre = 1; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 1;
Minv = @(v)feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,numcycles);

%See help Vcycle

%Parameters of GMRES iteration
tol   = 1e-8;
maxit = 50;

%GMRES iteration without preconditioner (too slow for large wavenumbers)
%tic
%[~,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
%time1 = toc;
 
% GMRES iteration with left SL preconditioner
tic
profile on
[x2,flag2,relres2,iter2,resvec2] = gmres(A,b,[],tol,maxit,Minv);
profile off
time2 = toc;

%GMRES iteration with right SL preconditioner
AMinv = @(v) A*feval(Minv,v);
tic
[x3,flag3,relres3,iter3,resvec3] = gmres(AMinv,b,[],tol,maxit);
x3 = feval(Minv,x3);
time3=toc;

%semilogy(1:(iter1(2)+1),resvec1'/resvec1(1),'b-+');
%hold on
semilogy(1:(iter2(2)+1),resvec2'/resvec2(1),'r-+');
hold on
semilogy(1:(iter3(2)+1),resvec3'/resvec3(1),'k-*');

%relative error with respect to exact solution
%relerr  = norm(x2-u_ex)/norm(u_ex)
%relerr2 = norm(x3-u_ex)/norm(u_ex)

%relerr = norm(A\b-u_ex)/norm(u_ex);

%%print results
%fprintf('Number of GMRES iterations (no preconditioning) %i\n',iter1(2));
fprintf('Number of GMRES iterations (+MG left  V-cycle preconditioning) %i\n',iter2(2));
fprintf('Number of GMRES iterations (+MG right V-cycle preconditioning) %i\n',iter3(2));


