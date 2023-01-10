function [x_sol] = Vcycle(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles)
%% VCYCLE Solves the Helmholtz/Poisson equation using a multigrid V-cycle.
%   Use: Vcycle_som(mg_mat,mg_split,restrict,interp,x0,b,npre,npos,w,smo,numcycles)
%   Input:
%       Output from the function mg_setup: 
%       - mg_mat:      cell array with multigrid matrices on all
%                      levels
%       - mg_split:    cell array with matrix splittings on all levels
%                      (for the smoothers)
%       - restrict:    cell array with restriction operators
%       - interp:      cell array with interpolation operators
%
%       x0:              initial guess
%       b:               right-hand side
%       smo:             smoother ('gs' for Gauss-Seidel, 'wjac' for w-Jacobi,'rbgs' for red-black Gauss Seidel in 2D)
%       npre, npos:      number of pre, post smoothing steps
%       w:               parameter for Jacobi iteration (set w=1 when using Gauss-Seidel)
%       numcycles:       number of V-cycles
%
%   Output       
%       x_sol:     solution computed with the multigrid V-cycle   
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 2.0, Sep 2016
%
% Notes:  The option 'liv' (livshits) is there to implement
%         Kaczmarcz relaxation, but this is not yet finished
%
%

%%    
    if length(mg_mat) == 1 %If on coarsest level, solve exactly
        %size(b)
        x_sol = mg_mat{1,1}\b;
        return;
        
    else
        x_sol = x0;
        
        for i=1:numcycles            
            %Presmoothing and computation of the residual
            %fprintf('Presmoothing with matrix of size %d\n',length(grid_matrices{1,1}));
            %Special case: Kaczmarcz relaxation
            if (length(mg_split{1})>4 && strcmp(smo,'liv')) 
                %Gauss-Seidel applied to normal equations
                rhs   = mg_mat{1}'*b;
                x_sol = smoother(mg_split{1}.Uk, mg_split{1}.Lk,...
                       mg_split{1}.Dk, mg_split{1}.P,rhs,x_sol,w,4,'wjac');                   
            elseif strcmp(smo,'liv')
                x_sol = smoother(mg_split{1}.U, mg_split{1}.L,...
                        mg_split{1}.D, mg_split{1}.P,b,x_sol,w,npre,'wjac');            
            else
               x_sol = smoother(mg_split{1}.U, mg_split{1}.L,...
                        mg_split{1}.D, mg_split{1}.P,b,x_sol,w,npre,smo);
            end
            
            res   = b-mg_mat{1}*x_sol;

            %Restriction of the residual to coarse grid
            fc   = restrict{1}*res; 
            nc   = length(fc); vc = zeros(nc,1);
            levs = length(mg_mat);
        
            %Calling Vcycle to solve the error equation
            if(levs >2)
                vc  = Vcycle(mg_mat(2:levs),mg_split(2:levs),restrict(2:levs-1),interp(2:levs-1),vc,fc,npre,npos,w,smo,1);
                %size(fc)
                
            elseif(levs==2)
                vc = Vcycle(mg_mat(2:2),mg_split(2:2),restrict(1:1),interp(1:1),vc,fc,npre,npos,w,smo,1);
                %size(fc)
            end
        
            %Prolongation and correction on current grid
            x_sol = x_sol + interp{1}*vc;
        
            %Postsmoothing  
            if (length(mg_split{1})>4 && strcmp(smo,'liv')) 
                %4 iter of Gauss-Seidel applied to normal equations
                rhs = mg_mat{1}'*b;
                x_sol = smoother(mg_split{1}.Uk, mg_split{1}.Lk,...
                       mg_split{1}.Dk,mg_split{1}.P,rhs,x_sol,w,4,'wjac');
                   
            elseif strcmp(smo,'liv')
                x_sol = smoother(mg_split{1}.U, mg_split{1}.L,...
                        mg_split{1}.D, mg_split{1}.P,b,x_sol,w,npos,'wjac');            
            else
               x_sol = smoother(mg_split{1}.U, mg_split{1}.L,...
                        mg_split{1}.D, mg_split{1}.P,b,x_sol,w,npos,smo);
            end
            
                         
        end
    end   

end

