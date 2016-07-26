%% Test multigrid Poisson 1D
npcc        = 5;                     %number of points in coarsest grid (1D) 
k           = 20; 
ppw         = 12; 
[npf,lev]   = npcc2npf(npcc,k,ppw);  %number of points in finest grid (1D)
bc          = 'dir';
b1 = 1; b2  = 0.5;
k  = 0;

%% Helmholtz and Shifted Laplacian Matrices and rhs
dim = 1;
A   = helmholtz(k,npf,bc);
b   = zeros(length(A),1); b(ceil(length(A)/4),1)=1;
z   = zeros(length(A),1); 

 
%% Multigrid Setup

% Construct matrices on all grids and interpolation operators
[grid_matrices,grid_smooth,restrict,interp] = mgsm_setup(A,lev,bc,dim);
 b  = zeros(length(A),1);
 x0 = sin(2*pi*randn(length(A),1));

% Parameters of V-cycle and Jacobi iteration
npre = 1; npos = 1; w = 2/3; smo = 'gs'; numcycles = 10;
r0   = norm(b-A*x0);

% Test of multigrid on 1D Poisson problem
% Use the profiler to evaluate performance of code
% profile on
% [xsol] = Vcycle(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
% profile off
% 
% r1  = b-A*xsol;
% resrat = norm(r1)/norm(r0);
% 
% r0=b-A*x0;
% %We compute the ratio between residuals 
% res_rat = zeros(numcycles-1,1);
% normres = zeros(numcycles,1); 
% err_rat = zeros(numcycles-1,1);

for i=1:numcycles
[xsol] = Vcycle(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,1);
r1     = b-A*xsol;
normres(i,1) = norm(r1);
res_rat(i,1) = norm(r1)/norm(r0);
err_rat(i,1) = norm(xsol)/norm(x0);
x0     = xsol; r0 = b-A*x0;
end


%% Test of Preconditioned GMRES

%Setting the SL preconditioner Minv
%Minv = @(v)feval(@Vcycle,SLgrid_matrices,restrict,interp,z,v,nu,mu,w);
% [L,U] = lu(M);
% Minv  = @(v) U\(L\v);

%Parameters of GMRES iteration
% tol   = 1e-8;
% maxit = length(A);
 
%GMRES iteration without preconditioner
% tic
% [~,flag1,relres1,iter1,resvec1] = ...
%   gmres(A,b,[],tol,maxit,[]);
% time1 = toc;

%GMRES iteration with SL preconditioner
% tic
% [x,flag2,relres2,iter2,resvec2] = ...
%    gmres(A,b,[],tol,maxit,Minv);
% time3 = toc;
%  
%semilogy(1:(max(iter2)+1),resvec2/resvec2(1),'r-+',1:(max(iter1)+1),resvec1/resvec1(1),'b-+');
