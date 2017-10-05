%Multigrid test of shifted Laplacian
clear global; 

%Parameters of Helmholtz equation and shifted Laplacian
dim      = 1;
k        = 50;
eps      = 0.5*k^2;   %Imaginary part of shift (for shifted Laplacian)
ppw      = 12;        %number of points per wavelength (fine grid)
npc      = 3; %number of points in the coarsest grid
op_type  = 'gal';
bc       = 'som';

%Number of points according to discretization rules
[npf,numlev] = fem_npc_to_npf(npc,k,ppw);

%Manual number of points
%numlev  = 10;
%[(2^(numlev-1))*npc,numlev];

h = 1/npf; grid = h*(1:1:npf); 
M = mass(npf,'mix');               %Mass matrix (for the norm)
norm2 = @(x)sqrt(abs(x'*M*x)); 

%%   Test of V-cycle on Helmholtz problem
% We test the multigrid solver on the Helmholtz problem
% Result: converges for approx k<=12, then diverges)
 
% Matrix hierarchy and right hand side 
f = ones(npf,1); f(npf) = 0.5; h = 1/npf; f = h*f; %(constant function f=1)
eps = 0;
[mg_mat,mg_split,restrict,interp] = mg_setupfem(k,eps,op_type,npc,numlev,dim);
u_exact = exact_sol(k,grid);
A = mg_mat{1}; 

% Parameters of V-cycle and smoother
 npre  = 2; npos = 1; w  = 2/3; numit = 10; smo = 'wjac';
 u0    = zeros(size(f));
 r0    = norm2(f-A*u0);
 u_ex  = A\f;
 res   = zeros(numit,1); 
 rat   = zeros(numit,1);
 
  for i=1:numit
      u_sol    = Vcyclefem(mg_mat,mg_split,restrict,...
                        interp,u0,f,npre,npos,w,smo,1);
      u0       = u_sol;
      res(i+1) = norm2(f-A*u0);
      rat(i)   = res(i+1)/res(i);
  end
  
 
%% Test of V-Cycle on shifted Laplace problem
% We test the multigrid solver on the shifted Laplace problem

% Matrix hierarchy and right hand side 
f   = ones(npf,1); f(npf)=0.5; h=1/npf; f=h*f; %(constant function f=1, note the scaling by h)
eps = 0.5*k^2;
[mg_mat_sl,mg_split_sl,restrict_sl,interp_sl] = mg_setupfem(k,eps,op_type,npc,numlev,dim);
M = mg_mat_sl{1};
u0 = zeros(length(M),1);

% Parameters of V-cycle and smoother
 npre = 2; npos = 1; w  = 2/3; numit = 10; smo = 'wjac';
 r0   = norm2(f-M*u0);
 u_ex = M\f;
 res_sl  = zeros(numit,1); 
 rat_sl  = zeros(numit,1);

 for i=1:numit
     u_sol    = Vcyclefem(mg_mat_sl,mg_split_sl,restrict,...
                       interp,u0,f,npre,npos,w,smo,1);
     u0       = u_sol;
     res_sl(i+1) = norm2(M*u_sol-f);
     rat_sl(i)   = res_sl(i+1)/res_sl(i);
 end

 
semilogy(res_sl)
title('Residual vs number of multigrid iterations');

rat
