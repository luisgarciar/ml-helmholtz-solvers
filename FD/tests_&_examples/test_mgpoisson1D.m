%% Test multigrid Poisson 1D
npc = 5;    %number of points in coarsest grid (1D) 
bc = 'dir';
m  = 3;

%wavenumber and imaginary shift of shifted Laplacian
k=10;  eps = 0; %Helmholtz problem
ppw = 12;   %number of points per wavelength

[npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)

k = 0;  eps = 0; %Poisson problem

%Exact solution (test problem)
u = @(x)  sin(m*pi*x);
f = @(x) ((m*pi)^2-k^2)*sin(m*pi*x);
h = 1/(npf+1);  grid = h*(1:1:npf)';
b = f(grid);    u_ex = u(grid);

%% Helmholtz and Shifted Laplacian Matrices and rhs
dim = 1;
A   = helmholtz(k,eps,npf,bc);
 
%% Multigrid Setup
% Construct matrices on all grids and interpolation operators
op_type = 'gal';
[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);
x0 = zeros(length(A),1);

% Parameters of V-cycle and Jacobi iteration
npre = 2; npos = 1; w = 2/3; smo = 'gs'; numcycles = 10;
r0   = norm(b-A*x0);

% Test of multigrid on 1D Poisson problem
% Use the profiler to evaluate performance of code
% profile on
 [xsol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
% profile off
% 
r1  = b-A*xsol;
relres = norm(r1)/norm(r0);
% 
for i=1:numcycles
[xsol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,1);
r1     = b-A*xsol;
normres(i,1) = norm(r1);
res_rat(i,1) = norm(r1)/norm(r0);
err_rat(i,1) = norm(xsol-x0)/norm(x0);
x0  = xsol; r0 = b-A*x0;
end

res_rat