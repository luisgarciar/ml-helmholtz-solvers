%% Solving Helmholtz problems with GMRES preconditioned by the Shifted Laplacian 
% 2-D Example, Dirichlet boundary conditions

clc
clear global;
close all;
npc = 1;    %number of interior points in coarsest grid in one dim 
bc  ='dir'; 
dim = 2;    %boundary conditions, dimension

%wavenumber and imaginary shift of shifted Laplacian
k   = 100;  eps = 0.5*k^2; %Helmholtz problem
ppw = 12;                 %number of points per wavelength
[npf,numlev] = fd_npc_to_npf(npc,k,ppw);  %number of points in finest grid (1D)

%% %Exact solution of Dirichlet problem
% u = @(x,y)  x.*(x-1).^2.*(y.^2).*(y-1);
% f = @(x,y) -k^2*x.*(y.^2).*((x-1).^2).*(y-1)...
%             -2*x.*((x-1).^2).*(3*y-1)...
%             -2*(y.^2).*(3*x-2).*(y-1);
% 
% %Helmholtz and Shifted Laplacian Matrices and rhs
% dim = 2;
% M   = helmholtz2(k,eps,npf,npf,bc);
% 
% npx = npf; npy=npf; hx = 1/(npx+1); hy = 1/(npy+1);
% xgrid = hx*(1:1:npx); ygrid = hy*(1:1:npy);
% np = npx*npy;
% [X,Y]= meshgrid(xgrid,ygrid);
% 
% u_ex = u(X,Y); u_ex = u_ex'; u_ex = reshape(u_ex,[np,1]); %exact solution
% b    = f(X,Y); b = b'; b = reshape(b,[np,1]);  %right hand side

%% Multigrid Setup
%npf = 63;    numlev=6;
profile on
A       = helmholtz2(k,0,npf,npf,bc);
op_type = 'gal'; %type of coarse operators (galerkin or rediscretized)
[mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npc,numlev,bc,dim);

% Setting the SL preconditioner Minv
% Parameters of V-cycle and Jacobi iteration
b    = ones(length(A),1);
x0   = zeros(size(b));
npre = 1; npos = 1; w = 0.5; smo = 'wjac'; numcycles = 1;
Minv_mg = @(v)feval(@Vcycle,mg_mat,mg_split,restrict,interp,x0,v,npre,npos,w,smo,numcycles);
AMinv_mg = @(v) A*feval(Minv_mg,v);


[L,U]   = lu(mg_mat{1});
Minv_ex  = @(v) U\(L\v);
AMinv_ex = @(v) A*feval(Minv_ex,v);


% %Parameters of GMRES iteration
tol   = 1e-6;
maxit = 400;

%GMRES iteration without preconditioner (too slow for large wavenumbers)
%tic
%[~,flag1,relres1,iter1,resvec1] = gmres(A,b,[],tol,maxit,[]);
%time1 = toc;
 
% GMRES iteration with left SL preconditioner
tic
[x2,flag2,relres2,iter2,resvec2] = gmres(AMinv_mg,b,[],tol,maxit);
profile off
time2 = toc;

%GMRES iteration with right SL preconditioner
AMinv = @(v) A*feval(Minv_mg,v);
tic
[x3,flag3,relres3,iter3,resvec3] = gmres(AMinv_ex,b,[],tol,maxit);
%x3 = feval(Minv,x3);
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
fprintf('Number of GMRES iterations (+MG right  V-cycle preconditioning) %i\n',iter2(2));
fprintf('Number of GMRES iterations (+MG exact preconditioning) %i\n',iter3(2));
