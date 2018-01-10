%% Test multigrid and GMRES, Sommerfeld BCs vs Dirichlet BCs
clear global; 
close all;
clc;

%% Parameters
k     = 60;       %wavenumber
ppw   = 15;       %min points per wavelength%
npcg  = 1;        %number of points in coarsest grid
dim   = 2;        %dimension
eps   = 0.5*k^2 ; %imaginary shift of shifted Laplacian

% Parameters of multigrid solver
npre = 1; npos = 1; w = 1/3; smo = 'liv';

%% Sommerfeld problem
bc           = 'som';      %boundary conditions
[npf,numlev] = fd_npc_to_npf(npcg,k,ppw); 
op_type      = 'gal';

%Multigrid setup for shifted Laplacian
[mg_mat_som,mg_split_som,restrict_som,interp_som] = mg_setup(k,eps,op_type,npcg,numlev,bc,dim);

%%
% right hand side and initial guess
M_som        = mg_mat_som{1}; 
%ex_sol_dir = ones(length(M_som),1);
%b_som      = M_som*ex_sol_dir;
b_som       = zeros(length(M_som),1); 
m           = floor(length(M_som)/2);
b_som(m,1)  = 1;
x0_som      = zeros(length(M_som),1); 

numcycles = 12;
% multigrid cycle for shifted Laplacian Sommerfeld
relres_som = zeros(numcycles+1,1); relres_som(1)=norm(b_som-M_som*x0_som); 
% relerr_som = zeros(numcycles+1,1); relerr_som(1)=norm(ex_sol_dir);

for i=1:numcycles
    [x_sol] = Vcycle(mg_mat_som,mg_split_som,restrict_som,interp_som,x0_som,b_som,npre,npos,w,smo,1);
    relres_som(i+1)=norm(b_som-M_som*x_sol);
    x0_som =x_sol;
end
 
relres_som = relres_som/relres_som(1);
factor_som = relres_som(2:length(relres_som))./relres_som(1:length(relres_som)-1);


% %% Dirichlet problem
% bc = 'dir';
% 
% % Multigrid setup for shifted Laplacian
% [mg_mat_dir,mg_split_dir,restrict_dir,interp_dir] = mg_setup(k,eps,op_type,npcg,numlev,bc,dim);
% 
% %Matrix, right hand side and initial guess
% M_dir      =  mg_mat_dir{1}; 
% %ex_sol_dir =  ones(length(M_dir),1);
% %b_dir      =  M_dir*ex_sol_dir;
% b_dir = zeros(length(M_dir),1); 
% m = floor(length(M_dir)/2);
% b_dir(m,1) = 1;
% x0  =  zeros(length(M_dir),1); 
% 
% %Setup of residuals and error
% relres_dir = zeros(numcycles+1,1); relres_dir(1) = norm(b_dir-M_dir*x0); 
% %relerr_dir = zeros(numcycles+1,1); relerr_dir(1) = norm(ex_sol_dir);
% 
% tic
% for i=1:numcycles
%     [x_sol] = Vcycle(mg_mat_dir,mg_split_dir,restrict_dir,interp_dir,x0,b_dir,npre,npos,w,smo,1);
%     relres_dir(i+1)=norm(b_dir-M_dir*x_sol);
%     x0 =x_sol;
% end
% time_dir = toc;  
% 
% relres_dir = relres_dir/relres_dir(1);
% factor_dir = relres_dir(2:length(relres_dir))./relres_dir(1:length(relres_dir)-1);
% 
% relres_som
% relres_dir
% 
% factor_som
% factor_dir


%% Test of MG + GMRES (Shifted Laplacian + Helmholtz)

%% Dirichlet problem
A_dir      = helmholtz2(k,0,npf,npf,bc);
%ex_sol_dir = ones(length(A_dir),1);
%b_dir      = A_dir*ex_sol_dir;
b_dir = zeros(length(M_dir),1); 
m = floor(length(M_dir)/2);
b_dir(m,1) = 1;
x0         = zeros(length(A_dir),1);
numcycles  = 2;
Minv_dir   = @(v)feval(@Vcycle,mg_mat_dir,mg_split_dir,restrict_dir,interp_dir,x0,v,npre,npos,w,smo,numcycles);
%[L_dir, U_dir] = lu(M_dir);
%Minv_dir = @(v) U_dir\(L_dir\v);
AMinv_dir = @(v)A_dir*feval(Minv_dir,v);


%GMRes parameters
tol   = 1e-8;
maxit = min(300,length(A_dir));

% profile on
% tic
% [x_lgdir,flag_lgdir,relres_lgdir,iter_lgdir,resvec_lgdir] = gmres(A_dir,b_dir,[],tol,maxit,Minv_dir);
% time_lgdir = toc;
% 
% tic
% [x_rgdir,flag_rgdir,relres_rgdir,iter_rgdir,resvec_rgdir] = gmres(AMinv_dir,b_dir,[],tol,maxit);
% time_rgdir = toc;
% 
% profile off

%% Sommerfeld Problem
bc = 'som';
A_som      = helmholtz2(k,0,npf,npf,bc);
%ex_sol_som = ones(length(A_som),1);
%b_som      = A_som*ex_sol_som;

b_som = zeros(length(A_som),1); b_som(floor(length(A_som)/2),1)=1;
x0 = zeros(length(A_som),1); 
numcycles = 1;
Minv_som = @(v) feval(@Vcycle,mg_mat_som,mg_split_som,restrict_som,interp_som,x0,v,npre,npos,w,smo,numcycles);
%[L_som, U_som] = lu(M_som);
%Minv_som = @(v) U_som\(L_som\v);
AMinv_som = @(v) A_som*feval(Minv_som,v);

%GMRes parameters
tol   = 1e-6;
maxit = min(200,size(A_som,1));

tic
[x_lgsom,flag_lgsom,relres_lgsom,iter_lgsom,resvec_lgsom] = gmres(A_som,b_som,[],tol,maxit,Minv_som);
time_lgsom = toc;

%Right-Preconditioned GMRES
tic
[x_rgsom,flag_rgsom,relres_rgsom,iter_rgsom,resvec_rgsom] = gmres(AMinv_som,b_som,[],tol,maxit);
time_rgsom = toc;


%% Plots of GMRES results
figure(1)
semilogy(0:(length(resvec_rgsom)-1), resvec_rgsom/resvec_rgsom(1), 'k-');
%hold on
%semilogy(0:(length(resvec_rgdir)-1), resvec_rgdir/resvec_rgdir(1), 'b-');
ylabel('relative residual')
xlabel('iteration')

legend('Sommerfeld BCs', 'Dirichlet BCs')
title(['1D Helmholtz with right CSL-preconditioner (k=',num2str(k),')'])

iter_rgsom
%iter_rgdir

time_rgsom
%time_rgdir

%% Plots of GMRES results
% figure(2)
% semilogy(0:iter_lgsom(2), resvec_lgsom/resvec_lgsom(1), 'k-');
% hold on
% semilogy(0:iter_lgdir(2), resvec_lgdir/resvec_lgdir(1), 'b-');
% ylabel('relative residual')
% xlabel('iteration')
% 
% legend('Sommerfeld BCs', 'Dirichlet BCs')
% title(['2D Helmholtz with left CSL-preconditioner (k=',num2str(k),')'])
% 
% iter_lgsom
% iter_lgdir
% 
% resvec_lgsom(end)/resvec_lgsom(1)
% resvec_lgdir(end)/resvec_lgdir(1)
% 
% time_lgsom
% time_lgdir

legend('Sommerfeld BCs', 'Dirichlet BCs')
title(['2D Helmholtz with left CSL-preconditioner (k=',num2str(k),')'])

iter_lgsom
iter_lgdir

resvec_lgsom(end)/resvec_lgsom(1)
resvec_lgdir(end)/resvec_lgdir(1)

time_lgsom
time_lgdir

