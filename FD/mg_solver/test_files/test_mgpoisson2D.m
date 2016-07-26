%Test multigrid Poisson
clear all;
npcc       = 5;                     %number of points in coarsest grid (1D) 
ppw        = 12;
k          = 20;
dim        = 2;
[npf,lev]  = npcc2npf(npcc,k,ppw);  %number of interior points in finest grid (1D)
bc         = 'dir';
k          = 0;
eps        = 0.5*k^2;


%% Poisson matrix and right hand side
npf = 100; 
npx = npf; npy = npf;  
%npx = 100; npy = 120; 
n=npx*npy;

A  = helmholtz2(k,0,npx,npy,bc);
%M  = shift_laplace2(k,eps,npx,npy,bc);

% For the exact solution u=sin(pi*2x)*sin(3*pi*y) the right hand side is
%   f = (-k^2+13*pi^2)*sin(2*pi*x)*sin(3*pi*y)
hx    = 1/(npx+1); gridx = hx*(1:1:npx); 
hy    = 1/(npy+1); gridy = hy*(1:1:npy);
[x,y] = meshgrid(gridx,gridy);
f     = @(x,y) (-k^2+13*pi^2)*(sin(2*pi*x).*sin(3*pi*y));
b     = feval(f,x,y); b = b'; b = reshape(b,[n,1]); %right hand side vector
u_ex  = (sin(2*pi*x).*sin(3*pi*y)); %exact solution
%figure(1); surf(x,y,u_ex);  %u_ex=u_ex';


% u_ex = reshape(u_ex,[n,1]); %exact solution
% u_d = A\b; 
% norm(u_ex-u_d)

%u_d=reshape(u_d',[npx,npy]); %discrete solution; 
%figure(2); surf(x,y,u_d);


% 
% 
% 
% %% Multigrid Setup
%  tic
%  [SLgrid_matrices,restrict,interp] = mgsm_setup(A,lev,bc,dim);
%  setuptime = toc;
% 
% %%
% %Multigrid Parameters
% nu = 2; mu = 2; w = 2/3;  
% 
% % Construct matrices on all grids and interpolation operators
% [grid_matrices,grid_smooth,restrict,interp] = mgsm_setup(A,lev,bc,dim);
%  x0 = sin(2*pi*randn(length(A),1));
% 
% % Parameters of V-cycle and Jacobi iteration
% npre = 2; npos = 2; w = 2/3; smo = 'gs'; numcycles = 50;
% r0   = norm(b-A*x0);
% 
% % Test of multigrid on 2D Poisson problem
% 
% % Use the profiler to evaluate performance of code
% profile on
% [xsol] = Vcycle(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,numcycles);
% profile off
% 
% r1  = b-A*xsol;
% resrat = norm(r1)/norm(r0);
% 
% 
% 
% % r0=b-A*x0;
% % %We compute the ratio between residuals 
% % res_rat = zeros(numcycles-1,1);
% % normres = zeros(numcycles,1); 
% % err_rat = zeros(numcycles-1,1);
% % 
% for i=1:numcycles
%     x0     = xsol; r0 = b-A*x0;
%     [xsol] = Vcycle(grid_matrices,grid_smooth,restrict,interp,x0,b,npre,npos,w,smo,1);
%     r1     = b-A*xsol;
%     normres(i,1) = norm(r1);
%     res_rat(i,1) = norm(r1)/norm(r0);
%     err_rat(i,1) = norm(xsol)/norm(x0);
% end
% 
% 
% %% Test of Preconditioned GMRES
% 
% %Setting the SL preconditioner Minv
% %Minv = @(v)feval(@Vcycle,SLgrid_matrices,restrict,interp,z,v,nu,mu,w);
% % [L,U] = lu(M);
% % Minv  = @(v) U\(L\v);
% 
% %Parameters of GMRES iteration
% % tol   = 1e-8;
% % maxit = length(A);
%  
% %GMRES iteration without preconditioner
% % tic
% % [~,flag1,relres1,iter1,resvec1] = ...
% %   gmres(A,b,[],tol,maxit,[]);
% % time1 = toc;
% 
% 
% %GMRES iteration with SL preconditioner
% % tic
% % [x,flag2,relres2,iter2,resvec2] = ...
% %    gmres(A,b,[],tol,maxit,Minv);
% % time3 = toc;
% %  
% %semilogy(1:(max(iter2)+1),resvec2/resvec2(1),'r-+',1:(max(iter1)+1),resvec1/resvec1(1),'b-+');
