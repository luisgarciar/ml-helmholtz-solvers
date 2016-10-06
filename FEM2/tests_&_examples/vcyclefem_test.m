%Multigrid test of shifted Laplacian
clear all; 

%Parameters of Helmholtz equation and shifted Laplacian
dim  = 1;
k    = 10;
eps  = 0.5*k^2;   %Imaginary part of shift (for shifted Laplacian)
ppw  = 20;        %number of points per wavelength (fine grid)
npc  = 4;         %number of points in the coarsest grid

[npf,numlev] = npc_to_npf(npc,k,ppw);
assert(k/npf<0.5,'grid is too coarse');

h = 1/npf; grid = h*(1:1:npf); 
A = helmholtzfem(k,npf,0);   %Helmholtz matrix
S = helmholtzfem(k,npf,eps); %Shifted Laplace matrix
M = mass(npf);               %Mass matrix (for the norm)
norm2 = @(x)sqrt(abs(x'*M*x)); 

%%   Test of V-cycle on Helmholtz problem
% % We test the multigrid solver on the Helmholtz problem
% % Result: converges for approx k<=12, then diverges)
% 
% % Matrix hierarchy and right hand side 
% f = ones(npf,1); f(npf)=0.5; h=1/npf; f=h*f; %(constant function f=1)
% [galerkin_matrices,galerkin_split,restrict,interp] = mg_setupfem(A,numlev,dim);
% u0 = zeros(length(A),1);
% u_exact = exact_sol(k,grid);
% 
% % Parameters of V-cycle and smoother
%  npre  = 2; npos = 2; w  = 2/3; numit = 10; smo = 'wjac';
%  r0    = norm2(f-A*u0);
%  u_ex  = A\f;
%  err   = zeros(numit,1); 
%  rat   = zeros(numit,1);
% 
%  for i=1:numit
%      u_sol    = Vcyclefem(galerkin_matrices,galerkin_split,restrict,...
%                        interp,u0,f,npre,npos,w,smo,1);
%      u0       = u_sol;
%      err(i+1) = norm2(u_sol-u_ex);
%      rat(i)   = err(i+1)/err(i);
%  end
%  
 
%% Test of V-Cycle on shifted Laplace problem
% We test the multigrid solver on the shifted Laplace problem

% Matrix hierarchy and right hand side 
f=ones(npf,1); f(npf)=0.5; h=1/npf; f=h*f; %(constant function f=1, note the scaling by h)
[galerkin_matrices,galerkin_split,restrict,interp] = mg_setupfem(S,numlev,dim);
u0 = zeros(length(S),1);

% Parameters of V-cycle and smoother
 npre = 2; npos = 2; w  = 1/3; numit = 10; smo = 'wjac';
 r0   = norm2(f-S*u0);
 u_ex = S\f;
 err  = zeros(numit,1); 
 rat  = zeros(numit,1);

 for i=1:numit
     u_sol    = Vcyclefem(galerkin_matrices,galerkin_split,restrict,...
                       interp,u0,f,npre,npos,w,smo,1);
     u0       = u_sol;
     err(i+1) = norm2(u_sol-u_ex);
     rat(i)   = err(i+1)/err(i);
 end

semilogy(err)
title('Error vs number of multigrid iterations');

rat
