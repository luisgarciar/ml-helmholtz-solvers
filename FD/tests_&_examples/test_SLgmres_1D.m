%% Test preconditioned GMRES for Helmholtz problems
npc = 3;    %number of interior points in coarsest grid in one dim 
bc = 'dir';

%wavenumber and imaginary shift of shifted Laplacian
k   = 100;  eps = 0.5*k^2; %Helmholtz problem
ppw = 12;   %number of points per wavelength
[npf,lev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)

%Exact solution (test problem)
% m = 3;
% u = @(x)  sin(m*pi*x);
% f = @(x) ((m*pi)^2-k^2)*sin(m*pi*x);
% h = 1/(npf+1);  grid = h*(1:1:npf)';
% b = f(grid);    u_ex = u(grid);






%% Helmholtz and Shifted Laplacian Matrices and rhs
dim = 1;
A  = helmholtz(k,0,npf,bc);  %set eps=0 for Helmholtz problem
M  = helmholtz(k,eps,npf,bc); %shifted Laplacian
b  = zeros(length(A),1); b(round(length(A)/2),1)=1;
 
%% Multigrid Setup
% Construct matrices on all grids and interpolation operators
[grid_matrices,grid_smooth,restrict,interp] = mg_setup(M,lev,bc,dim);
x0 = zeros(length(M),1);
 

%% Test of Preconditioned GMRES
%Setting the MG preconditioner Minv
%Parameters of V-cycle and Jacobi iteration
npre = 2; npos = 2; w = 2/3; smo = 'wjac'; numcycles = 1;
Minv = @(v)feval(@Vcycle,grid_matrices,grid_smooth,restrict,interp,x0,v,npre,npos,w,smo,1);

%Parameters of GMRES iteration
tol   = 1e-8;
maxit = length(A);
 
%GMRES iteration without preconditioner
 tic
 [x1,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
 time1 = toc;

% GMRES iteration with left SL preconditioner
tic
[x2,flag2,relres2,iter2,resvec2] = gmres(A,b,[],tol,maxit,Minv);
time2 = toc;

% GMRES iteration with right SL preconditioner
AMinv = @(v) A*feval(Minv,v);
tic
[x3,flag3,relres3,iter3,resvec3] = gmres(AMinv,b,[],tol,maxit);
time3=toc;

iter1
iter2
iter3

%Note: Comparing the timings for this problem does not tell much
%Because the solution is real and using the shifted Laplacian
%requires complex arithmetic
time1
time2
time3

semilogy(1:(max(iter2)+1),resvec2/resvec2(1),'r-+',...
         1:(max(iter1)+1),resvec1/resvec1(1),'b-+',...
         1:(max(iter3)+1),resvec3/resvec3(1),'k-*')