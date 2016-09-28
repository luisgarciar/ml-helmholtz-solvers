function [galerkin_matrices,galerkin_split,restrict,interp] = mg_setup_som(A,numlev,bc,dim)
%% MG_SETUP_SOM: Constructs a hierarchy of Galerkin coarse grid matrices,
%            splitting of the Galerkin matrices, and interpolation operators for a Helmholtz/Laplace PDE problem   
%  Use:
% [galerkin_matrices,galerkin_split,restrict,interp] = mg_setup_som(A,numlev,bc,dim)
%
%  Input: 
%  A:         Matrix on finest grid
%  numlev:    Total number of levels (number of coarse grids)
%  bc:        Type of boundary conditions:      
%             'dir' for homogeneous dirichlet bc's
%             'som' for sommerfeld bc's 
%  dim:        Dimension (1 or 2)
%
%  Output:
%  grid_matrices: Cell array with Galerkin matrices 
%                 (grid_matrices{i}: Galerkin matrix at level i)   
%
%  grid_split: Cell array with splitting of Galerkin matrices
%              (upper, lower and diagonal part) to be applied in
%              smoothing steps 
%              (grid_split{i}.U, grid_split{i}.L, grid_split{i}.D)
% 
%  restrict:   Cell array with restriction operators
%              (restrict_op{i}:restriction operator from i to i+1)
%
%  interp:     Cell array with interpolation operators
%              (interp_op{i}: interpolation operator from i+1 to i)
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%
%  Version 2.0, Sep 2016
%  Works with Sommerfeld boundary conditions
%
%%  Construction of restriction, interpolation and Galerkin matrices
s  = length(A);
[npc,npf] = size2npc_som(s,dim,bc);

restrict          = cell(numlev-1,1);    %restrict{i}: fw restriction grid i to grid i+1 
interp            = cell(numlev-1,1);    %interp{i}: lin interp grid i+1 to grid i
galerkin_matrices = cell(numlev,1);      %grid_matrices{i}: Galerkin matrix at level i
galerkin_split    = cell(numlev,1);      %grid_split{i}: matrix splittings needed for smoothers i

%Level 1 
i=1;
galerkin_matrices{i} = A;
restrict{i}         = fwrestriction_som(npf,dim,bc); %fw restriction, grid 1 to grid 2 
interp{i}           = lininterpol_som(npc,dim,bc);      %lin interp, grid 2 to grid 1
galerkin_split{i}.U = sparse(triu(A,1));         %matrix splitting of A
galerkin_split{i}.L = sparse(tril(A,-1));
galerkin_split{i}.D = spdiags(diag(A),0,length(A),length(A));
galerkin_split{i}.P = perm_rb_som(length(galerkin_matrices{i}));

% [s1,s2] = size(grid_matrices{i});
% [r1,r2] = size(restrict{i});
% [i1,i2] = size(interp{i});

% formatText1 = 'Level %d: Number of fine points %d, number of coarse points %d\n';
% formatText2 = 'Level %d: Size of linear system %dx%d, size of restriction matrix %dx%d, size of prolongation matrix %dx%d\n';
% 
% sprintf(formatText1,i,npf,npc)
% sprintf(formatText2,i,s1,s2,r1,r2,i1,i2)
              
% if dim==2
%     grid_smooth{1}.P=rb_reorder(npf);
% end

if numlev==1
    return;
end


for i=2:numlev-1
    npf = npc;  npc = round((npf-1)/2);  
    galerkin_matrices{i} = sparse(restrict{i-1}*galerkin_matrices{i-1}*interp{i-1}); %Coarse grid Galerkin matrix
    galerkin_split{i}.U=sparse(triu(galerkin_matrices{i},1));
    galerkin_split{i}.L=sparse(tril(galerkin_matrices{i},-1));
    galerkin_split{i}.D=spdiags(diag(galerkin_matrices{i}),0,length(galerkin_matrices{i}),length(galerkin_matrices{i}));
    galerkin_split{i}.P=perm_rb_som(length(galerkin_matrices{i}));

    if i<numlev
        restrict{i}  = fwrestriction_som(npf,dim,bc); %fw restriction, grid i to grid i+1 
        interp{i}    = lininterpol_som(npc,dim,bc);   %lin interp, grid i+1 to grid i
    end

%             if dim==2
%                 grid_smooth{i}.P=sparse(rb_reorder(npf));
%             end

end

npf = npc;  
galerkin_matrices{numlev} = sparse(restrict{numlev-1}*galerkin_matrices{numlev-1}*interp{numlev-1});
galerkin_split{numlev}.U  = sparse(triu(galerkin_matrices{numlev},1));
galerkin_split{numlev}.L  = sparse(tril(galerkin_matrices{numlev},-1));
galerkin_split{numlev}.D  = spdiags(diag(galerkin_matrices{numlev}),0,length(galerkin_matrices{numlev}),length(galerkin_matrices{numlev}));
galerkin_split{numlev}.P  = perm_rb_som(length(galerkin_matrices{numlev}));

%       if dim==2 %Red black permutation matrix
%           grid_smooth{numlev}.P=rb_reorder(npf);
%       end
end       


