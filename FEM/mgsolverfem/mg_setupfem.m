function [mg_mat,mg_split,restrict,interp] = mg_setupfem(k,eps,op_type,npcc,numlev,dim)
%% MGSM_SETUPFEM: Constructs a hierarchy of Galerkin coarse grid matrices,
% smoother splittings and interpolation operators for a Helmholtz problem
% discretized with P1 finite elements
%
%  Use:
% [galerkin_matrices,galerkin_split,restrict,interp] = mg_setupfem(A,numlev,dim)
%
%  Input: 
%  k:
%  eps:     Parameters of Helmholtz problem (k=0 for Poisson) and
%           shifted Laplacian 
%           (see 'help helmholtz2' for info on choosing eps)
%
%  op_type:   Type of coarse grid operators
%             'rd'  for rediscretized operators on coarse levels
%             'gal' for Galerkin operators on coarse levels
%
%  npcc:      number of interior points on coarsest grid in 1D
%             (x=0 and interior points)
%
%  numlev:    Total number of levels (number of coarse grids + 1)
%
%  Output:
%  mg_mat: Cell array with Galerkin matrices 
%                     (grid_matrices{i}: Galerkin matrix level i)   
%
%  mg_split: Cell array with upper, lower and diagonal part of
%                  Galerkin matrices to be applied in smoothing steps
%                  (galerkin_split{i}.U, galerkin_split{i}.L, galerkin_split{i}.D)
% 
%  restrict:       Cell array with restriction operators
%                  (restrict_op{i}:restriction operator from i to i+1)
%
%  interp:      Cell array with interpolation operators
%              (interp_op{i}:restriction operator from i+1 to i)
%
% Reference: F. Ihlenburg and I. Babuska, Finite element solution of the 
% Helmholtz equation with high wave number Part I: The h-version of the FEM
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 2.0, Oct 2016
%
% To Do : Add 2D case

%%  Construction of restriction and interpolation matrices

%% to be fixed
npf = npc_numlev_to_npf_fem(npcc,numlev); %add this function%
npc = round(npf/2);  %check this for fem discretizations!!
%%

restrict = cell(numlev-1,1);    %restrict{i}: fw restriction grid i to grid i+1 
interp = cell(numlev-1,1);    %interp{i}: lin interp grid i+1 to grid i
mg_mat = cell(numlev,1);      %grid_matrices{i}: Galerkin matrix at level i
mg_split = cell(numlev,1);      %grid_smoothers{i}: matrix splittings needed for smoothers i

%Level 1 
mg_mat{1} = helmholtzfem(k,eps,npf,bc);
restrict{1} = fwrestrictionfem(npf,dim);  %fw restriction, grid 1 to grid 2 
interp{1}   = lininterpolfem(npc,dim);    %lin interp, grid 2 to grid 1

mg_split{1}.U = sparse(triu(mg_mat{1},1));  %matrix splitting of mg_mat{}
mg_split{1}.L = sparse(tril(mg_mat{1},-1));
mg_split{1}.D = spdiags(diag(mg_mat{1}),0,length(mg_mat{1}),length(mg_mat{1}));
mg_split{1}.P = speye(length(mg_mat{mg_mat{1}})); %red-black permutation matrix 

if numlev==1
    return;
end

npf = npc;   npc = round(npf/2);         

for i=2:numlev-1
    
    if i<numlev      
        restrict{i} = fwrestrictionfem(npf,dim); %fw restriction, grid i to grid i+1
        interp{i}   = lininterpolfem(npc,dim);   %lin interp, grid i+1 to grid i        
    end
    
    switch op_type
        case 'rd'  %Coarse matrix by rediscretization
           switch dim
               case 1
                   mg_mat{i} = helmholtzfem(k,eps,npf,bc);
               case 2
                   %mg_mat{i} = helmholtz2(k,eps,npf,npf,bc); 
                   error ('invalid dimension') %add FEM 2D later%
               otherwise
                   error('invalid dimension')
           end
           
        case 'gal' %Galerkin coarse matrix
           mg_mat{i} = sparse(restrict{i-1}*mg_mat{i-1}*interp{i-1});

        otherwise
        error('Invalid operator type');      
    end    
    
    mg_split{i}.U  = sparse(triu(mg_mat{i},1));
    mg_split{i}.L  = sparse(tril(mg_mat{i},-1));
    mg_split{i}.D  = spdiags(diag(mg_mat{i}),0,length(mg_mat{i}),length(mg_mat{i}));
    mg_split{i}.P  = speye(length(mg_mat{i}));    
    npf = npc;   npc = round(npf/2);            

end
            
%Coarsest Level
switch op_type
        case 'rd'  %Coarse matrix by rediscretization
            switch dim
               case 1                  
                   mg_mat{numlev} = helmholtzfem(k,eps,npf,bc);  
               case 2
                   error('invalid dimension') %fix later
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
mg_split{numlev}.P =  speye(length(mg_mat{numlev})); %fix later

end       


