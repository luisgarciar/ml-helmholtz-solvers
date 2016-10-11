%% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian 
%  1-D Example, Dirichlet boundary conditions
%  Discretized with finite differences

npc = 1;    %number of interior points in coarsest grid in one dim 
bc = 'dir';
%wavenumber and imaginary shift of shifted Laplacian
k   = 50;  eps = 0.5*k^2; %Helmholtz problem
ppw = 20;   %number of points per wavelength
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)

%Exact solution (test problem, not good for GMRES, eigenfunction!)
% m = 3;
% u = @(x)  sin(m*pi*x);
% f = @(x) ((m*pi)^2-k^2)*sin(m*pi*x);
% h = 1/(npf+1);  grid = h*(1:1:npf)';
% b = f(grid);    u_ex = u(grid);

%Exact solution
 u = @(x) x.*(x-1).^2; 
 f = @(x) -6.*x+4-k^2*x.*(x-1).^2;
 h = 1/(npf+1);  grid = h*(1:1:npf)';
 b = f(grid);    u_ex = u(grid);

%% Helmholtz and Shifted Laplacian Matrices and rhs
dim = 1;
A  = helmholtz(k,0,npf,bc);  %set eps=0 for Helmholtz problem
%M  = helmholtz(k,eps,npf,bc); %shifted Laplacian
 
%% Multigrid Setup
% Construct matrices on all grids and interpolation operators
op_type = 'gal';
[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);

x0 = zeros(length(A),1);
 
%% Test of Preconditioned GMRES
%Setting the MG preconditioner Minv
%Parameters of V-cycle and Jacobi iteration
npre = 2; npos = 1; w = 2/3; smo = 'wjac'; numcycles = 1;
Minv = @(v)feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,1);

%Parameters of GMRES iteration
tol   = 1e-10;
maxit = length(A);
 
%GMRES iteration without preconditioner
 tic
 [x1,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
 time1 = toc;

%GMRES iteration with left SL preconditioner
tic
[x2,flag2,relres2,iter2,resvec2] = gmres(A,b,[],tol,maxit,Minv);
time2 = toc;

%GMRES iteration with right SL preconditioner
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

semilogy(1:(iter2(2)+1),resvec2'/resvec2(1),'r-+')
hold on
semilogy(1:(iter1(2)),resvec1'/resvec1(1),'b-+')
semilogy(1:(iter3(2)+1),resvec3'/resvec3(1),'k-*');

relerr1=norm(x1-u_ex)/norm(u_ex)
relerr2=norm(x2-u_ex)/norm(u_ex)
%relerr3=norm(x3-M*u_ex)/norm(M*u_ex)