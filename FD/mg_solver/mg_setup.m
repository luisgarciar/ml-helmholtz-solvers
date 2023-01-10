function [mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npcc,numlev,bc,dim)
%% MG_SETUP: Constructs a hierarchy of coarse grid operators,
%            splittings and intergrid operators for a Helmholtz/Laplace PDE problem
%            in 2D
%  Use:
%
% [mg_mat,mg_split,restrict,interp] = mg_setup(k,eps,op_type,npcc,numlev,bc,dim)
%
%  Input: 
%
%  k,eps:     Parameters of Helmholtz problem (k=0 for Poisson) and
%             shifted Laplacian 
%             (see 'help helmholtz2' for info on choosing eps)
%
%  op_type:   Type of coarse grid operators
%             'rd'  for rediscretized operators on coarse levels
%             'gal' for Galerkin operators on coarse levels
%
%  npcc:      number of interior points on coarsest grid in 1D
%
%  numlev:    Total number of levels (number of coarse grids + 1)
%
%  bc:        Type of boundary conditions:      
%             'dir' for homogeneous dirichlet bc's
%             'som' for sommerfeld bc's
%             'som1' for sommerfeld bc's with 1st order discretization  
%
%  dim:        Dimension (1 or 2)
%
%  Output:
%
%  mg_mat:     Cell array with  matrices at all levels
%              (mg_mat{i}: matrix at level i)   
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

%[npc,npf] = size2npc(s,dim,bc);
%s  = length(A);
npf = npc_numlev_to_npf(npcc,numlev);
restrict  = cell(numlev-1,1);    %restrict{i}: fw restriction grid i to grid i+1 
interp    = cell(numlev-1,1);    %interp{i}: lin interp grid i+1 to grid i
mg_mat    = cell(numlev,1);      %grid_matrices{i}: Galerkin matrix at level i
mg_split  = cell(numlev,1);      %grid_split{i}: matrix splittings needed for smoothers i

%For Kaczmarcz
hmax = 1/(npcc+1);
hlevs = hmax*fliplr(2.^-(0:1:(numlev-1)));

%kczlevel is the level at which Kaczmarcz relaxation should be performed
[~,kczlevel] = min(abs(k*hlevs - 1.25));  
if length(kczlevel)>1
    kczlevel=kczlevel(1);
end

%Level 1 
switch dim
    case 1
    mg_mat{1} = helmholtz(k,eps,npf,bc);
    case 2
    %npf
    mg_mat{1} = helmholtz2(k,eps,npf,npf,bc);
    otherwise
        error('invalid dimension')
end

npc = round((npf-1)/2);

restrict{1}   = fwrestriction(npf,dim,bc);  %fw restriction, grid 1 to grid 2 
interp{1}     = lininterpol(npc,dim,bc);    %lin interp, grid 2 to grid 1
mg_split{1}.U = sparse(triu(mg_mat{1},1));  %matrix splitting of A
mg_split{1}.L = sparse(tril(mg_mat{1},-1));
mg_split{1}.D = spdiags(diag(mg_mat{1}),0,length(mg_mat{1}),length(mg_mat{1}));
mg_split{1}.P = speye(length(mg_mat{1}));

if dim==2
    mg_split{1}.P = perm_rb(length(mg_mat{1}));
end

if numlev==1
    return;
end

npf = npc; npc= round((npf-1)/2);

for i=2:numlev-1   
    if i<numlev      
        restrict{i} = fwrestriction(npf,dim,bc); %fw restriction, grid i to grid i+1
        interp{i}   = lininterpol(npc,dim,bc);   %lin interp, grid i+1 to grid i        
    end
    
    switch op_type
        case 'rd'  %Coarse matrix by rediscretization
           switch dim
               case 1
                   mg_mat{i} = helmholtz(k,eps,npf,bc);
               case 2
                   mg_mat{i} = helmholtz2(k,eps,npf,npf,bc); 
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
    
    %If on Kaczmarcz level, save the splittings for
    %Kaczmarcz relaxation also
    if i == kczlevel
        mg_split{i}.Uk = sparse(triu(mg_mat{i}'*mg_mat{i},1));
        mg_split{i}.Lk = sparse(tril(mg_mat{i}'*mg_mat{i},1));
        mg_split{i}.Dk = sparse(diag(mg_mat{i}'*mg_mat{i},1));
    end
    
    switch dim
        case 1
            mg_split{i}.P = eye(npf); 
        case 2
            mg_split{i}.P = eye(npf); 
    end    
    npf = npc; npc= round((npf-1)/2);
  
end

%Coarsest Level
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


