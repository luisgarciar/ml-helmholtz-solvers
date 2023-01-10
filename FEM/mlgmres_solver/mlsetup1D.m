function [ml_mat,ml_prec,ml_prec_split,restrict,interp] = mlsetup1D(k1,eps1,k2,eps2,op_type,npcc,numlev,dim,bc)
% MLSETUP1D: Constructs a hierarchy of Galerkin coarse grid matrices,
% smoother splittings and interpolation operators for a 2D Helmholtz/shifted
% Laplace problem discretized with P1 finite elements on a uniform square
% grid 
%
%  Use:
% [ml_mat,ml_prec,ml_prec_split,restrict,interp] = mlsetup1D(k1,eps1,k2,eps2,op_type,npcc,numlev,dim,bc)
%
%  Input: 
%  k1, eps1:     Parameters of Helmholtz problem  
%           (see 'help helmholtz2' for info on choosing eps)
%  k2, eps2:     Parameters of Shifted Laplace preconditioner
%
%  op_type:   Type of coarse grid operators
%             'rd'  for rediscretized operators on coarse levels
%             'gal' for Galerkin operators on coarse levels
%
%  npcc:      number of interior points on coarsest grid in 1D          
%  numlev:    Total number of levels (number of coarse grids + 1)
%  pde:       Structure with the data of the pde problem
%  option:    structure with options
%            - option.twolevel: If 'true', the function only returns finest
%               and coarsest grid
%
%  Output:
%  ml_mat: Cell array with Galerkin matrices 
%          (grid_matrices{i}: Galerkin matrix level i)   
%
%  ml_prec: Cell array with preconditioners at level i
%           
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 1.0 Sep 2019
%%

[ml_mat,~,restrict,interp]  = mg_setupfem_1D(k1,eps1,op_type,npcc,numlev,dim,bc);
[ml_prec,ml_prec_split,~,~] = mg_setupfem_1D(k2,eps2,op_type,npcc,numlev,dim,bc);

end       


