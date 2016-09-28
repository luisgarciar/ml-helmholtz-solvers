function Z = lininterpol(npc, dim)
%% LININTERPOL Constructs the matrix corresponding to the linear interpolation
%  operator from a coarse grid with npc interior points.
%  
%   Use:    Z = lininterpol(npc,dim)  
%
%   Input: 
%       npc:    number of 1D interior points in coarse grid
%       dim:    dimension (1 or 2)
%       
%   Output:
%       Z:      interpolation matrix
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
%          Version 1.0, Jun 2016
%          Corrected and checked 2D
%
%   To Do: Add Other Boundary conditions
%
%  References: 
%  A Multigrid Tutorial, Briggs, Henson, McCormick, 2000, Chap 3
%  Computational Science and Engineering, Strang, 2007, Chap 7 
%
%%
switch dim
    case 1
        npf = 2*(npc)+1; %length of fine grid vectors
        z   = zeros(npf,1); z(1:3,1) = [1;2;1];
        Z   = gallery('circul',z)';
        Z   = 0.5*sparse(Z(:,1:2:(npf-2)));
        
    case 2
        npf = 2*(npc)+1; %length of fine grid vectors (1D)
        z   = zeros(npf,1); z(1:3,1) = [1;2;1];
        Z   = gallery('circul',z)';
        Z   = 0.5*sparse(Z(:,1:2:(npf-2))); %1D operator
        Z   = kron(Z,Z);  %2D operator        
end
        
        
end

