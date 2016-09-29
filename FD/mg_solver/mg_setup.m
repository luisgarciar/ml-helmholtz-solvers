function [mg_mat,mg_split,restrict,interp] = mg_setup(A,k,eps,op_type,numlev,bc,dim)
%% MG_SETUP: Constructs a hierarchy of coarse grid operators,
%            splittings and intergrid operators for a Helmholtz/Laplace PDE problem   
%  Use:
%
% [mg_mat,mg_split,restrict,interp] = mg_setup(A,k,eps,op_type,numlev,bc,dim)
%
%  Input: 
%  A:         Matrix on finest grid
%
%  k,eps:     Parameters of Helmholtz problem (k=0 for Poisson) and
%             shifted Laplacian 
%             (see 'help helmholtz2' for info on choosing eps)
%
%  op_type:   Type of coarse grid operators
%             'rd'  for rediscretized operators on coarse levels
%             'gal' for Galerkin operators on coarse levels
%
%  numlev:    Total number of levels (number of coarse grids + 1)
%
%  bc:        Type of boundary conditions:      
%             'dir' for homogeneous dirichlet bc's
%             'som' for sommerfeld bc's 
%
%  dim:        Dimension (1 or 2)
%
%  Output:
%
%  mg_mat:     Cell array with multigrid matrices 
%                 (mg_mat{i}: matrix at level i)   
%
%  mg_split:   Cell array with splitting of mg_mat matrices
%              (upper, lower and diagonal part) to be applied in
%              smoothing steps 
%              (mg_split{i}.U, mg_split{i}.L, mg_split{i}.D)
%               mg_split{i}.P: Red-Black permutation matrix
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
%  Added option for rediscretized operators on coarse grids
%
%%  Construction of restriction, interpolation and Galerkin matrices
s  = length(A);
[npc,npf] = size2npc(s,dim,bc);

restrict  = cell(numlev-1,1);    %restrict{i}: fw restriction grid i to grid i+1 
interp    = cell(numlev-1,1);    %interp{i}: lin interp grid i+1 to grid i
mg_mat    = cell(numlev,1);      %grid_matrices{i}: Galerkin matrix at level i
mg_split  = cell(numlev,1);      %grid_split{i}: matrix splittings needed for smoothers i

%Level 1 
mg_mat{1}     = A;
restrict{1}   = fwrestriction(npf,dim,bc); %fw restriction, grid 1 to grid 2 
interp{1}     = lininterpol(npc,dim,bc);   %lin interp, grid 2 to grid 1
mg_split{1}.U = sparse(triu(A,1));         %matrix splitting of A
mg_split{1}.L = sparse(tril(A,-1));
mg_split{1}.D = spdiags(diag(A),0,length(A),length(A));
mg_split{1}.P = eye(length(A));

if dim==2
    mg_split{1}.P = perm_rb(length(mg_mat{1}));
end

if numlev==1
    return;
end

for i=2:numlev-1
    npf = npc;  npc = round((npf-1)/2);     
    switch op_type
        case 'rd'  %Coarse matrix by rediscretization
           switch dim
               case 1
                   mg_mat{i} = helmholtz(k,eps,npf,bc);
                   size(mg_mat{i});
               case 2
                   mg_mat{i} = helmholtz2(k,eps,npf,npf,bc); 
                   size(mg_mat{i});
               otherwise
                   error('invalid dimension')
           end
           
        case 'gal' %Galerkin coarse matrix
           mg_mat{i} = sparse(restrict{i-1}*mg_mat{i-1}*interp{i-1}); 
        otherwise
        error('Invalid operator type');      
    end    
    
    mg_split{i}.U = sparse(triu(mg_mat{i},1));
    mg_split{i}.L = sparse(tril(mg_mat{i},-1));
    mg_split{i}.D = spdiags(diag(mg_mat{i}),0,length(mg_mat{i}),length(mg_mat{i}));
    
    switch dim
        case 1
            mg_split{i}.P = eye(npf); 
        case 2
            mg_split{i}.P = perm_rb(length(mg_mat{i}));
    end
            
    if i<numlev
        restrict{i} = fwrestriction(npf,dim,bc); %fw restriction, grid i to grid i+1 
        interp{i}   = lininterpol(npc,dim,bc);   %lin interp, grid i+1 to grid i
    end

end

%Coarsest Level
npf = npc;  

switch op_type
        case 'rd'  %Coarse matrix by rediscretization
            switch dim
               case 2
                   mg_mat{numlev} = helmholtz2(k,eps,npf,npf,bc);  
               case 1
                   mg_mat{numlev} = helmholtz(k,eps,npf,bc);
                otherwise
                   error('Invalid dimension')
            end
        case 'gal' %Galerkin coarse matrix
            mg_mat{numlev} = sparse(restrict{numlev-1}*mg_mat{numlev-1}*interp{numlev-1});
        otherwise
        error('Invalid operator type');
end
    
mg_split{numlev}.U  = sparse(triu(mg_mat{numlev},1));
mg_split{numlev}.L  = sparse(tril(mg_mat{numlev},-1));
mg_split{numlev}.D  = spdiags(diag(mg_mat{numlev}),0,length(mg_mat{numlev}),length(mg_mat{numlev}));

switch dim
    case 1
        mg_split{numlev}.P = eye(npf); 
    case 2
        mg_split{numlev}.P = perm_rb(length(mg_mat{numlev}));
end

end       


