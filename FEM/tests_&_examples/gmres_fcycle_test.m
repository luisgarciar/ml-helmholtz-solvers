%Test of GMRES preconditioned by the shifted Laplacian for 
%the 1-D Helmholtz problem
%
%       -u''- k^2*u = f in (0,1)
%        u(0)       = 0   
%        u'(1)-ku(1)= 0

clear all
close all
clc

k    = 60;        %wavenumber
eps  = 0.5*k^2;   %imaginary shift of shifted Laplacian
npc  = 2;         %number of points in the coarsest grid
ppw  = 40;        %number of points per wavelength
dim  = 1;

[npf,numlev] = npc_to_npf(npc,k,ppw);  %number of fine grid points and levels
h = 1/npf; grid = h*(1:1:npf)'; 

A = helmholtzfem(k,npf,0);   %Helmholtz matrix
S = helmholtzfem(k,npf,eps); %Shifted Laplace matrix
M = mass(npf);               %Mass matrix (for the norm)
norm2 = @(x)sqrt(abs(x'*M*x)); 

%right hand side %(constant function f=1, note the scaling by h)
f = ones(npf,1); f(npf)=0.5; h=1/npf; f=h*f; 
u_ex   = exact_sol(k,grid); %analytic solution
u_dir  = A\f;    %solution computed by direct method (to check the discr. error)
relerr_dir = norm2(u_ex-u_dir)/norm2(u_ex); %Relative error in the L2 norm

%Including the endpoint x=0 in the grid and the solutions
%grid=[0;grid]; u_ex=[0;u_ex]; u_dir=[0;u_dir];

%plot of discrete vs analytical solution
figure(1)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_dir),'b-');
title('Analytic solution vs Discrete Solution');

%Solving Au=f with GMRES and no preconditioner
restart = []; %% number of iter before restart, [] = no restart.
tol     = 1e-8;   %% tolerance for relative residual
maxit   = npf;

tic 
[u_gmres,~,~,iter_gmres,resvec_gmres] = gmres(A, f, restart, tol, npf) ;
time_gmres = toc;

%Relative error of GMRES solution
relerr_gmres = norm2(u_gmres-u_ex)/norm2(u_ex);

%Setting up the multigrid preconditioner for the shifted Laplacian
[galerkin_matrices,galerkin_split,restrict,interp] = mg_setupfem(S,numlev,dim);
u0 = zeros(length(S),1);

npre = 2; npos = 2; w = 2/3; numit = 1; smo = 'wjac';
Sinv = @(v) Fcyclefem(galerkin_matrices,galerkin_split,restrict,...
            interp,v,npre,npos,w,smo);
Sinv2 = @(v) Vcyclefem(galerkin_matrices,galerkin_split,restrict,...
            interp,u0,v,npre,npos,w,smo,numit);
                 
ASinv = @(v) A*Sinv(v); %Right Preconditioned matrix
         
%If we use right preconditioning we solve in two steps: ASinv*v=f, u=Sinv*v;
profile clear      
tic   
profile on
[u_rpgmres,~,~,iter_rpgmres,resvec_rpgmres] = gmres(ASinv, f, restart, tol, npf);
profile off
time_rpgmres= toc;

u_rpgmres = S\u_rpgmres; %we only compute u to check the error (note the amplification in the error!)
relerr_rpgmres = norm2(u_rpgmres-u_ex)/norm2(u_ex);

%Alternative: Left preconditioning
%One can also use left preconditioning and pass the preconditioner as
%an argument to the gmres function (left prec. is default in MATLAB gmres)
tic   
[u_lfpgmres,~,~,iter_lfpgmres,resvec_lfpgmres] = gmres(A, f, restart, tol,npf,Sinv);
time_lfpgmres= toc;
relerr_lfpgmres = norm2(u_lfpgmres-u_ex)/norm2(u_ex);

tic   
[u_lvpgmres,~,~,iter_lvpgmres,resvec_lvpgmres] = gmres(A, f, restart, tol,npf,Sinv2);
time_lvpgmres= toc;
relerr_lvpgmres = norm2(u_lvpgmres-u_ex)/norm2(u_ex);

figure(2)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_rpgmres),'b-'); 
title('Analytic solution vs Iterative Solution, Right preconditioning');
legend('Analytic solution','GMRES solution (right prec.)')


figure(3)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_lfpgmres),'b-'); 
title('Analytic solution vs Iterative Solution, Left preconditioning');
legend('Analytic solution','GMRES solution (left prec.)')

figure(4)
plot(grid,real(u_ex),'r-');
hold on
plot(grid,real(u_lvpgmres),'b-'); 
title('Analytic solution vs Iterative Solution, Left preconditioning - Vcycle');
legend('Analytic solution','GMRES solution (left prec.)')

figure(5)
semilogy(1:length(resvec_gmres),resvec_gmres/resvec_gmres(1),'r',...
         1:length(resvec_rpgmres),resvec_rpgmres/resvec_rpgmres(1),'b',...
         1:length(resvec_lfpgmres),resvec_lfpgmres/resvec_lfpgmres(1),'k',...
         1:length(resvec_lvpgmres),resvec_lvpgmres/resvec_lvpgmres(1),'g');
          
legend('No Preconditioning','Right Preconditioning','Left F-Preconditioning','Left V-Preconditioning');
title('Relative residuals of GMRES iteration');


fprintf('Relative error direct vs analytic solution %d\n',relerr_dir);
fprintf('Relative error iterative (right prec) vs analytic solution %d\n',relerr_rpgmres);
fprintf('Relative error iterative (left  prec-Fcycle) vs analytic solution %d\n\n',relerr_lfpgmres);
fprintf('Relative error iterative (left  prec-Vcycle) vs analytic solution %d\n\n',relerr_lvpgmres);


fprintf('Number of GMRES iterations (no preconditioning) %i\n',iter_gmres(2));
fprintf('Number of GMRES iterations (+MG right preconditioning) %i\n',iter_rpgmres(2));
fprintf('Number of GMRES iterations (+MG left F-cycle preconditioning) %i\n\n',iter_lfpgmres(2));
fprintf('Number of GMRES iterations (+MG left V-cyle preconditioning) %i\n\n',iter_lvpgmres(2));




fprintf('Iteration time GMRES (no preconditioning)  %d\n',time_gmres);
fprintf('Iteration time GMRES (+MG right preconditioning) %d\n',time_rpgmres);
fprintf('Iteration time GMRES (+MG left F-cycle preconditioning) %d\n\n',time_lfpgmres);
fprintf('Iteration time GMRES (+MG left V-cycle preconditioning) %d\n\n',time_lvpgmres);

profile viewer