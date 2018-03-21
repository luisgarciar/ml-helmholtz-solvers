%Multigrid test of shifted Laplacian
clear global; 

%Parameters of Helmholtz equation and shifted Laplacian
k          = 80;
factoreps  = 0.5;
poweps     = 2;
eps        = factoreps*k^poweps;   %Imaginary part of shift (for shifted Laplacian)
ppw        = 0.5;                  %number of points per wavelength (fine grid)
npcc       = 4;                    %number of points in the coarsest grid
op_type    = 'gal';
bc         = 'som';

%Number of points according to discretization rules
[~,numlev] = fd_npc_to_npf(npcc,k,ppw);  %number of points in finest grid (1D)

pdeSL      = helmholtz2Dconstantwndata(k,factoreps,poweps);

%%   Test of V-cycle on shifted Laplace problem
%We test the multigrid solver on the shifted Laplace problem
 
%Matrix hierarchy and right hand side 
[mg_mat_fem,mg_split_fem,restr_fem,interp_fem]= mg_setupfem_2D(npcc,numlev,pdeSL);
Afem = mg_mat_fem{1};

% Parameters of V-cycle and smoother
 npre = 1; npos = 1; w  = 0.7; numit = 20; smo = 'wjac';
 u_ex = randn(length(Afem),1);
 f    = Afem*u_ex;
 u0   = sparse(length(Afem),1);
 r0   = norm(f-Afem*u0);
 res_fem  = zeros(numit,1); res_fem(1)=r0;
 rat  = zeros(numit,1);
 
 profile on
  for i=1:numit
      u_sol    = Wcycle(mg_mat_fem,mg_split_fem,restr_fem,...
                           interp_fem,u0,f,npre,npos,w,smo,1);
      u0       = u_sol;
      res_fem(i+1) = norm(f-Afem*u0);
      rat(i)   = res_fem(i+1)/res_fem(i);
  end
  profile off
  
  
  relres = res_fem/res_fem(1);
  relres

  
  
  
  
  