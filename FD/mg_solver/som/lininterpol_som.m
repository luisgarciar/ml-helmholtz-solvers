function Z = lininterpol_som(npc, dim,bc)
%% LININTERPOL_SOM Constructs the matrix corresponding to the (bi)linear 
%interpolation operator from a coarse grid with npc interior points.
%  
%   Use:    Z = lininterpol_som(npc,dim,bc)  
%
%   Input: 
%       npc:  number of 1D interior points in coarse grid
%       dim:  dimension (1 or 2)
%       bc:  'dir' for homogeneous Dirichlet boundary conditions
%            'som' for first order Sommerfeld boundary conditions 
%       
%   Output:
%       Z:  interpolation matrix, size npf x npc (1D Dirichlet)
%                                      npf^2 x npc^2 (2D Dirichlet)
%                                     (npf+2)^2 x (npc+2)^2 (2D Sommerfeld)
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 2.0, Jun 2016
%
%  Dirichlet and Sommerfeld boundary conditions
%
%  References: 
%  A Multigrid Tutorial, Briggs, Henson, McCormick, 2000, Chap 3
%  Computational Science and Engineering, Strang, 2007, Chap 7 
%  Multigrid, Trottenberg, Oosterlee, Sch?ller, 2001, Chap 2
%
%%
switch dim
    case 1
        npf = 2*(npc)+1; %length of fine grid vectors
        z   = zeros(npf,1); z(1:3,1) = [1;2;1];
        Z   = gallery('circul',z)';
        Z   = 0.5*sparse(Z(:,1:2:(npf-2)));
        
    case 2
        switch bc
            case 'dir'
                npf = 2*(npc)+1; %length of fine grid vectors (1D)
                z   = zeros(npf,1); z(1:3,1) = [1;2;1];
                Z   = gallery('circul',z)';
                Z   = 0.5*sparse(Z(:,1:2:(npf-2))); %1D operator
                Z   = kron(Z,Z);  %2D operator 
                
            case 'som'
               Z=4*fwrestriction_som(npf,dim,bc)';
        end
        
        
end

