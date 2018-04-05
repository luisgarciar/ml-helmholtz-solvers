%Multigrid test of shifted Laplacian
clear global; 

%Parameters of Helmholtz equation and shifted Laplacian
k          = 100;
factoreps  = 1;
poweps     = 2;
eps        = factoreps*k^poweps;    %Imaginary part of shift (for shifted Laplacian)
ppw        = 0.5;                   %number of points per wavelength (fine grid)
npcc       = 10;                    %number of points in the coarsest grid
op_type    = 'gal';
bc         = 'som';

%Number of points according to discretization rules
[~,numlev] = fd_npc_to_npf(npcc,k,ppw);  %number of points in finest grid (1D)
pdeSL      = helmholtz2Dconstantwndata(k,factoreps,poweps);


%%   Test of V-cycle on shifted Laplace problem
%We test the multigrid solver on the shifted Laplace problem
 
%Matrix hierarchy and right hand side 
[mg_mat,mg_split,restr,interp]= mg_setupfem_2D(npcc,numlev,pdeSL);
A = mg_mat{1};

% Parameters of V-cycle and smoother
 npre = 1; npos = 1; w  = 0.5; numit = 20; smo = 'wjac';
 u_ex = ones(length(A),1);
 f    = A*u_ex;
 u0   = sparse(length(A),1);
 r0   = norm(f-A*u0);
 res  = zeros(numit,1); 
 rat  = zeros(numit,1);
 
 res(1) = r0;

profile on
tic 
  for i=1:numit
      u_sol    = Vcycle(mg_mat,mg_split,restr,...
                           interp,u0,f,npre,npos,w,smo,1);
      u0       = u_sol;
      res(i+1) = norm(f-A*u0);
      rat(i)   = res(i+1)/res(i);
  end
  time_mg = toc
  profile off
  
  
  res/res(1)
  
