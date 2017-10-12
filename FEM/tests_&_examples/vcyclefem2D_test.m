%Multigrid test of shifted Laplacian
clear global; 

%Parameters of Helmholtz equation and shifted Laplacian
k          = 50;
factoreps  = 10;
poweps     = 2;
eps        = factoreps*k^poweps;   %Imaginary part of shift (for shifted Laplacian)
ppw        = 15;                   %number of points per wavelength (fine grid)
npcc       = 10;                    %number of points in the coarsest grid
op_type    = 'gal';
bc         = 'som';

%Number of points according to discretization rules
[~,numlev] = fd_npc_to_npf(npcc,k,ppw);  %number of points in finest grid (1D)
pdeSL   = helmholtz2Dconstantwndata(k,factoreps,poweps);


%%   Test of V-cycle on shifted Laplace problem
%We test the multigrid solver on the shifted Laplace problem
 
%Matrix hierarchy and right hand side 
[mg_mat,mg_split,restr,interp]= mg_setupfem_2D(npcc,numlev,pdeSL);
A = mg_mat{1};

% Parameters of V-cycle and smoother
 npre = 1; npos = 1; w  = 2/3; numit = 20; smo = 'wjac';
 u_ex = ones(length(A),1);
 f    = A*u_ex;
 u0   = sparse(length(A),1);
 r0   = norm(f-A*u0);
 res  = zeros(numit,1); 
 rat  = zeros(numit,1);
 
  for i=1:numit
      u_sol    = Vcyclefem(mg_mat,mg_split,restr,...
                        interp,u0,f,npre,npos,w,smo,1);
      u0       = u_sol;
      res(i+1) = norm(f-A*u0);
      rat(i)   = res(i+1)/res(i);
  end
  res
  rat
