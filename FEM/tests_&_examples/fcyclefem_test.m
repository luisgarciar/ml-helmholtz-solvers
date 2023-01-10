%Multigrid test of shifted Laplacian
clear all; close all; clc

%Parameters of Helmholtz equation and shifted Laplacian
dim  = 1;
k    = 100;
eps  = 0.5*k^2;   %Imaginary part of shift (for shifted Laplacian)
ppw  = 20;        %number of points per wavelength (fine grid)
npc  = 4;         %number of points in the coarsest grid

[npf,numlev] = npc_to_npf(npc,k,ppw);
assert(k/npf<0.5,'grid is too coarse');

k    = 20;
eps  = 0.5*k^2;  

h = 1/npf; grid = h*(1:1:npf); 
A = helmholtzfem(k,npf,0);   %Helmholtz matrix
S = helmholtzfem(k,npf,eps); %Shifted Laplace matrix
M = mass(npf);               %Mass matrix (for the norm)
norm2 = @(x)sqrt(abs(x'*M*x)); 


%% Test of F-Cycle on Helmholtz problem
% We test the multigrid solver on the Helmholtz problem (small k)

% Matrix hierarchy and right hand side 
f=ones(npf,1); f(npf)=0.5; h=1/npf; f=h*f; %(constant function f=1, note the scaling by h)
[galerkin_matrices,galerkin_split,restrict,interp] = mg_setupfem(S,numlev,dim);
u0 = zeros(length(A),1);
u_dir  = S\f;  %discrete solution with direct method
u_ex   = exact_sol(k,grid)'; %analytic solution
relerr_dir = norm2(u_ex-u_dir)/norm2(u_ex); %Relative error in the L2 norm

% Parameters of F-cycle and smoother
 npre   = 1; npos = 1; w  = 1/3; smo = 'gs'; numit = 100;
 r0     = norm2(f-S*u0);
 err    = zeros(numit,1); 
 relres = zeros(numit,1); 
 rat    = zeros(numit,1);
 relres(1)=1;
 
 for i=1:numit
     u_sol       = Fcyclefem(galerkin_matrices,galerkin_split,restrict,...
                       interp,f,npre,npos,w,smo);
     u0          = u_sol;
     relres(i+1) = norm2(f-S*u_sol);
     rat(i)      = relres(i+1)/relres(i);
 end

close all; 
figure(1)
semilogy(relres)
title('Relative residual vs number of multigrid iterations');

figure(2)
plot(grid,real(u_ex),'r-');
title('Analytic solution');
figure(3)
plot(grid,real(u_sol),'b-');
title('Multigrid solution');

rat
