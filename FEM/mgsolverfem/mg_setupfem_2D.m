function [mg_mat,mg_split,restr,interp] = mg_setupfem_2D(npcc,numlev,pde)
%% MG_SETUPFEM_2D: Constructs a hierarchy of Galerkin coarse grid matrices,
% smoother splittings and interpolation operators for a 2D Helmholtz/shifted
% Laplace problem discretized with P1 finite elements on a uniform square
% grid 
%
% Requires the iFEM package!! 
%
%  Use:
% [mg_mat,mg_split,restrict,interp] = mg_setupfem_2D(npcc,numlev,pde)
%
%
%  Input: 
%
%  npcc:      number of interior points on coarsest grid in 1D          
%  numlev:    Total number of levels (number of coarse grids + 1)
%  pde:       Structure with the data of the pde problem
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
%
%  interp:      Cell array with interpolation operators
%              (interp_op{i}:restriction operator from i+1 to i)
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

restr = flip(Res);    restr  = restr(1:numlev-1);
interp  = flip(Pro);  interp = interp(2:numlev);

mg_split = cell(numlev,1);    

for i=1:numlev
    
    %matrix splitting of mg_mat{i}
    mg_split{i}.U = sparse(triu(mg_mat{i},1));  
    mg_split{i}.L = sparse(tril(mg_mat{i},-1));
    mg_split{i}.D = spdiags(diag(mg_mat{i}),0,length(mg_mat{i}),length(mg_mat{i}));

end


end       


