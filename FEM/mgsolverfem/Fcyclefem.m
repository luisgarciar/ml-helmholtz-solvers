function [x_sol] = Fcyclefem(galerkin_matrices,galerkin_split,restrict_op,interp_op,rhs,npre,npos,w,smo)
%% FCYCLEFEM Solves the Helmholtz equation using a multigrid F-cycle.
%
%   Use: Fcyclefem(galerkin_matrices,galerkin_split,restrict_op,interp_op,x0,b,npre,npos,w,smo,numcycles)
%
%   Input:

%       Data from the function mgsmsetup: 
%         - galerkin_matrices:  cell array with matrices on all
%                               levels
%         - galerkin_split:     cell array with matrix splittings for smoothers
%         - restrict_op:        cell array with restriction operators
%         - interp_op:          cell array with interpolation operators
%
%       x0:              initial guess
%       rhs:             right-hand side
%       smo:             smoother ('gs' for Gauss-Seidel, 'wjac' for w-Jacobi,'rbgs' for red-black Gauss Seidel in 2D)
%       npre:            number of presmoothing steps
%       npos:            number of postsmoothing steps
%       w:               parameter for Jacobi iteration (set w=1 when using Gauss-Seidel)
%       numcycles:       number of V-cycles
%
%   Output       
%       x_sol:     solution computed with the multigrid F-cycle   
%
%   Reference: 
%   W. Briggs, V.E. Henson & S. McCormick, A Multigrid Tutorial, p.43
%
%   Author: Luis Garcia Ramos,          
%           Institut fur Mathematik, TU Berlin
%           Version 1.0, Jun 2016
%
%%  

    if length(galerkin_matrices) == 1 %If on coarse level, solve exactly
        x_sol = galerkin_matrices{1}\rhs;
        %x_sol = smoother(galerkin_split{1}.U,galerkin_split{1}.L,galerkin_split{1}.D,...
                        % galerkin_split{1}.P,rhs,zeros(size(rhs)),w,npre,smo);
        return;
        
    else %restrict to coarse grid and solve recursively
        rhsc  = restrict_op{1}*rhs; %rhs_coarse
        levs  = length(galerkin_matrices);
        xc    = Fcyclefem(galerkin_matrices(2:levs),galerkin_split(2:levs),restrict_op(2:levs-1),interp_op(2:levs-1),...
                          rhsc,npre,npos,w,smo);
    end
    
    x_sol = interp_op{1}*xc;
    x_sol = Vcyclefem(galerkin_matrices,galerkin_split,restrict_op,...,
                    interp_op,x_sol,rhs,npre,npos,w,smo,1);
                       

end

