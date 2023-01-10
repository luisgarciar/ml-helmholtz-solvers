%% Test of GMRES preconditioned by the shifted Laplacian for 
%the 1-D Helmholtz problem
%
%       -u''- k^2*u = f in (0,1)
%        u(0)       = 0   
%        u'(1)-ku(1)= 0

clear global; 
close all;
clc;

k    = 50;       %wavenumber
eps  = 0.5*k^2;   %imaginary shift of shifted Laplacian
npcg = 2;         %number of points in the coarsest grid
ppw  = 20;        %number of points per wavelength
dim  = 1;

[npf_fem,numlev_fem] = fem_npc_to_npf(npcg,k,ppw);  %number of fine grid points and levels
h = 1/npf_fem; grid  = h*(1:1:npf_fem)'; 

op_type = 'gal'; %galerkin coarse operators

A = helmholtzfem(k,npf_fem,0);   %Helmholtz matrix
S = helmholtzfem(k,npf_fem,eps); %Shifted Laplace matrix
M = mass(npf_fem);               %Mass matrix (for the norm)
norm2 = @(x)sqrt(abs(x'*M*x)); 

%right hand side %(constant function f=1, note the scaling by h)
f = ones(npf_fem,1); f(npf_fem)=0.5; h=1/npf_fem; rhs=h*f; 
u_ex   = exact_sol(k,grid); %analytic solution
u_dir  = A\rhs;    %solution computed by direct method (to check the discr. error)
relerr_dir = norm2(u_ex-u_dir)/norm2(u_ex); %Relative error in the L2 norm

%Including the endpoint x=0 in the grid and the solutions
%grid=[0;grid]; u_ex=[0;u_ex]; u_dir=[0;u_dir];

%plot of discrete vs analytical solution
figure(1)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_dir),'b-');
title('Analytic solution vs Discrete Solution');

%% 
%Solving Au=f with GMRES and no preconditioner
restart = []; %% number of iter before restart, [] = no restart.
tol     = 1e-8;   %% tolerance for relative residual
maxit   = npf_fem;

tic 
[u_gmres,~,~,iter_gmres,resvec_gmres] = gmres(A, rhs, restart, tol, npf_fem) ;
time_gmres = toc;

%Relative error of GMRES solution
relerr_gmres = norm2(u_gmres-u_ex)/norm2(u_ex);

%Setting up the multigrid preconditioner for the shifted Laplacian
[mg_mat,mg_split,restrict,interp] = mg_setupfem(k,eps,op_type,npcg,numlev_fem,dim);
u0 = zeros(length(S),1);

%size(mg_mat{1})

npre = 2; npos = 2; w = 2/3; numvcycles = 1; smo = 'wjac';
Sinv = @(v) Vcyclefem(mg_mat,mg_split,restrict,...
                      interp,u0,v,npre,npos,w,smo,numvcycles);
                               
ASinv = @(v) A*Sinv(v); %Right Preconditioned matrix
         
%If we use right preconditioning we solve in two steps: ASinv*v=f, u=Sinv*v;
profile clear      
tic   
profile on
[u_rpgmres,~,~,iter_rpgmres,resvec_rpgmres] = gmres(ASinv, rhs, restart, tol, npf_fem);
profile off
time_rpgmres= toc;

u_rpgmres = S\u_rpgmres; %we only compute u to check the error (note the amplification in the error!)
relerr_rpgmres = norm2(u_rpgmres-u_ex)/norm2(u_ex);

%Alternative: Left preconditioning
%One can also use left preconditioning and pass the preconditioner as
%an argument to the gmres function (left prec. is default in MATLAB gmres)
tic   
[u_lpgmres,~,~,iter_lpgmres,resvec_lpgmres] = gmres(A, rhs, restart, tol,npf_fem,Sinv);
time_lpgmres= toc;
relerr_lpgmres = norm2(u_lpgmres-u_ex)/norm2(u_ex);

figure(2)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_rpgmres),'b-'); 
title('Analytic solution vs Iterative Solution, Right preconditioning');
legend('Analytic solution','GMRES solution (right prec.)')

figure(3)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_lpgmres),'b-'); 
title('Analytic solution vs Iterative Solution, Left preconditioning');
legend('Analytic solution','GMRES solution (left prec.)')

figure(4)
semilogy(1:length(resvec_gmres),resvec_gmres/resvec_gmres(1),'r',...
         1:length(resvec_rpgmres),resvec_rpgmres/resvec_rpgmres(1),'b',...
          1:length(resvec_lpgmres),resvec_lpgmres/resvec_lpgmres(1),'k');
legend('No Preconditioning','Right Preconditioning','Left Preconditioning');
title('Relative residuals of GMRES iteration');


fprintf('Relative error direct vs analytic solution %d\n',relerr_dir);
fprintf('Relative error iterative (right prec) vs analytic solution %d\n',relerr_rpgmres);
fprintf('Relative error iterative (left  prec) vs analytic solution %d\n\n',relerr_lpgmres);

fprintf('Number of GMRES iterations (no preconditioning) %i\n',iter_gmres(2));
fprintf('Number of GMRES iterations (+MG right preconditioning) %i\n',iter_rpgmres(2));
fprintf('Number of GMRES iterations (+MG left preconditioning) %i\n\n',iter_lpgmres(2));

fprintf('Iteration time GMRES (no preconditioning)  %d\n',time_gmres);
fprintf('Iteration time GMRES (+MG right preconditioning) %d\n',time_rpgmres);
fprintf('Iteration time GMRES (+MG left preconditioning) %d\n\n',time_lpgmres);

profile viewer