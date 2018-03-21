function [mg_mat,mg_split,restrict,interp] = mg_setupfem_2D(npcc,numlev,pde,option)
%% MG_SETUPFEM_2D: Constructs a hierarchy of Galerkin coarse grid matrices,
% smoother splittings and interpolation operators for a 2D Helmholtz/shifted
% Laplace problem discretized with P1 finite elements on a uniform square
% grid 
%
% Requires the iFEM package!! 
%
%  Use:
% [mg_mat,mg_split,restr,interp] = mg_setupfem_2D(npcc,numlev,pde)
%
%  Input: 
%
%  npcc:      number of interior points on coarsest grid in 1D          
%  numlev:    Total number of levels (number of coarse grids + 1)
%  pde:       Structure with the data of the pde problem
%  option:    structure with options
%         - option.twolevel: If 'true', the function only returns finest
%           and coarsest grid
%
%  Output:
%  mg_mat: Cell array with Galerkin matrices 
%                     (grid_matrices{i}: Galerkin matrix level i)   
%
%  mg_split: Cell array with upper, lower and diagonal part of
%                  Galerkin matrices to be applied in smoothing steps
%                  (galerkin_split{i}.U, galerkin_split{i}.L, galerkin_split{i}.D)
% 
%  restrict:   Cell array with restriction operators
%              (restrict_op{i}:restriction operator from i to i+1)
%  If option.twolevel = 'true', restrict is a matrix of size according to 
%  finest and coarsest grids
%
%  interp:      Cell array with interpolation operators
%              (interp_op{i}:restriction operator from i+1 to i)
%  If option.twolevel = 'true', interp is a matrix of size according to 
%  finest and coarsest grids, and interp = restrict'
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 1.0 Sep 2017
%
%%

%npf = npc_numlev_to_npf_fem(npcc,numlev); %add this function%

%coarse mesh
h = 1/(npcc+1);
[node,elem] = squaremesh([0 1 0 1],h);

%refining the mesh numlev times
for j = 1:numlev-1
    [node,elem] = uniformrefine(node,elem);
end

%Find boundary nodes
[~,bdEdge,~] = findboundary(elem);
    
%Sets Sommerfeld boundary conditions on all boundary edges
bdFlag = setboundary(node,elem,'ABC');

%pde data
%pde    = helmholtz2Dconstantwndata(k,factoreps,poweps);
[eqn,~] = helmholtz2Dfem(node,elem,pde,bdFlag,bdEdge);

A = eqn.A;
n = size(A,1);
b = sparse(n,1);
option.solver = 'NO';

[~,~,Ai,~,~,Res,Pro,~] = mg(A,b,elem,option);
assert(length(Ai)==numlev, 'error: incorrect number of levels');

mg_mat = flip(Ai);

restrict = flip(Res);  restrict  = restrict(1:numlev-1);
interp   = flip(Pro);  interp    = interp(2:numlev);

mg_split = cell(numlev,1);    

for i=1:numlev
    
    %matrix splitting of mg_mat{i}
    mg_split{i}.U = sparse(triu(mg_mat{i},1));  
    mg_split{i}.L = sparse(tril(mg_mat{i},-1));
    mg_split{i}.D = spdiags(diag(mg_mat{i}),0,length(mg_mat{i}),length(mg_mat{i}));
    mg_split{i}.P = speye(length(mg_mat{i}));

end

%If only coarsest and finest levels are needed
if isfield(option,'twolevel') 
    if option.twolevel==true 
    two_lev{1} = mg_mat{1};
    two_lev{2} = mg_mat{end};
    mg_mat = two_lev;
    
    R = restrict{1};
    
    for i = 2:length(restrict)
        R = restrict{i}*R;
    end    
    
    restrict  = R;   
    interp    = restrict';
    
    [m,n] = size(restrict); %Restriction operator from fine to coarse grid
    assert(m==length(mg_mat{2}),'incorrect size of restriction operator');
    assert(n==length(mg_mat{1}),'incorrect size of restriction operator');
    end
    
end

end
       


