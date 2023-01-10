function [x_sol] = Vcyclefem(galerkin_matrices,galerkin_split,restrict_op,interp_op,x0,b,npre,npos,w,smo,numcycles)
%% VCYCLEFEM Solves the Helmholtz equation using a multigrid V-cycle.
%
%   Use: Vcyclefem(galerkin_matrices,galerkin_split,restrict_op,interp_op,x0,b,npre,npos,w,smo,numcycles)
%
%   Input:
%       Output from the function mg_setupfem or mg_setupfem_2D: 
%         - galerkin_matrices:  cell array with matrices on all
%                               levels
%         - galerkin_split:     cell array with matrix splittings for smoothers
%         - restrict_op:        cell array with restriction operators
%         - interp_op:          cell array with interpolation operators
%
%       x0:              Initial guess
%       b:               right-hand side
%       smo:             smoother ('gs' for Gauss-Seidel, 'wjac' for w-Jacobi,'rbgs' for red-black Gauss Seidel in 2D)
%       npre:            number of presmoothing steps
%       npos:            number of postsmoothing steps
%       w:               parameter for Jacobi iteration (set w=1 when using Gauss-Seidel)
%       numcycles:       number of V-cycles
%
%   Output       
%       x_sol:     solution computed with the multigrid V-cycle   
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 1.0, Jun 2016
%
%%    

if numcycles==1

    if length(galerkin_matrices) == 1 %If on coarse level, solve exactly
        x_sol = galerkin_matrices{1}\b;
        return;
    else
        x_sol = x0;
                    
            %Presmoothing and computation of the residual
            %fprintf('Presmoothing with matrix of size %d\n',length(grid_matrices{1,1}));
            x_sol = smoother(galerkin_split{1}.U, galerkin_split{1}.L,...
                             galerkin_split{1}.D, galerkin_split{1}.P,b,x_sol,w,npre,smo);
            res   = b-galerkin_matrices{1}*x_sol;

            %Restriction of the residual to coarse grid
            rc   = restrict_op{1}*res; 
            nc   = length(rc); vc = zeros(nc,1);
            levs = length(galerkin_matrices);
        
            %Calling Vcycle to solve the error equation
            if(levs >2)
                vc  = Vcyclefem(galerkin_matrices(2:levs),galerkin_split(2:levs),restrict_op(2:levs-1),interp_op(2:levs-1),vc,rc,npre,npos,w,smo,1);
                %size(fc)
                
            elseif(levs==2)
                vc = Vcyclefem(galerkin_matrices(2:2),galerkin_split(2:2),restrict_op(1:1),interp_op(1:1),vc,rc,npre,npos,w,smo,1);
                %size(fc)
            end
        
            %Prolongation and correction on current grid
            x_sol = x_sol + interp_op{1}*vc;
        
            %Postsmoothing
            x_sol = smoother(galerkin_split{1}.U,galerkin_split{1}.L,...
                             galerkin_split{1}.D,galerkin_split{1}.P,b,x_sol,w,npos,smo);           
    end 
    
else
    x_init = x0;
    for i=1:numcycles
        x_sol = Vcyclefem(galerkin_matrices,galerkin_split,restrict_op,interp_op,x_init,b,npre,npos,w,smo,1);
        x_init = x_sol;
    end
 
end

