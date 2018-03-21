function Z = lininterpol(npc,dim,bc)
%% LININTERPOL Constructs the matrix corresponding to the (bi)linear 
%  interpolation operator from a coarse grid with npc interior points.
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
%       Z:  interpolation matrix, size npf x npc     (1D Dirichlet)
%                                      npf^2 x npc^2 (2D Dirichlet)
%                                     (npf+2)^2 x (npc+2)^2 (2D Sommerfeld)
%
%  Author: Luis Garcia Ramos, 
%          Institut fur Mathematik, TU Berlin
% 
%Version 2.0, Sep 2016:
%Works with Dirichlet and Sommerfeld boundary conditions
%
%TO DO: Add option for different number of points in x and y direction in
%       2D
%
%  References: 
%  A Multigrid Tutorial, Briggs, Henson, McCormick, 2000, Chap 3
%  Computational Science and Engineering, Strang, 2007, Chap 7 
%  Multigrid, Trottenberg, Oosterlee, Schiller, 2001, Chap 2
%%
switch dim
    case 1
        switch bc
            case 'dir'
                npf = 2*(npc)+1; %length of fine grid vectors
                z   = zeros(npf,1); z(1:3,1) = [1;2;1];
                Z   = gallery('circul',z)';
                Z   = 0.5*sparse(Z(:,1:2:(npf-2)));
                
            case 'som'
                npcc = npc+2; npff = 2*npcc-1;                 
                u = 0.5*ones(npcc,1);
                odd = (1:2:npff);  %indices of coarse grid points
                even= (2:2:npff-1);
               
                %interpolation in 1D
                Z  = sparse(npff,npcc); 
                Z1 = spdiags([u u],[1 0],npcc-1,npcc);
                Z(odd,:) = speye(npcc); %points on both coarse & fine grids
                Z(even,:) = Z1;    %remaining points
                                
            otherwise
                error('invalid boundary conditions');
        end
        
    case 2
        npf = 2*(npc)+1; %number of interior points in fine grid (1D)
        
        switch bc
            case 'dir'                
                z   = zeros(npf,1); z(1:3,1) = [1;2;1];
                Z   = gallery('circul',z)';
                Z   = 0.5*sparse(Z(:,1:2:(npf-2))); %1D operator
                Z   = kron(Z,Z);  %2D operator 
                
            case 'som'               
               %modify this lines later to add option for
               %diff. number of points on x and y
               npcix  = npc; npciy = npc;
               npccx = npcix+2; npccy= npciy+2;
               npffx = 2*npccx-1; npffy = 2*npccy-1; 
               
               u = 0.5*ones(npccx,1);
               odd = (1:2:npffx);  %indices of coarse grid points
               even= (2:2:npffx-1);
               
               %interpolation in 1D
               Zx  = sparse(npffx,npccx); 
               Z1 = spdiags([u u],[1 0],npccx-1,npccx);
               Zx(odd,:) = speye(npccx); %coarse grid points
               Zx(even,:) = Z1;
                             
               Z = kron(Zx,Zx); %2D Operator
               
               case 'som1'               
               %modify this lines later to add option for
               %diff. number of points on x and y
               npcix  = npc; npciy = npc;
               npccx = npcix+2; npccy= npciy+2;
               npffx = 2*npccx-1; npffy = 2*npccy-1; 
               
               u = 0.5*ones(npccx,1);
               odd = (1:2:npffx);  %indices of coarse grid points
               even= (2:2:npffx-1);
               
               %interpolation in 1D
               Zx  = sparse(npffx,npccx); 
               Z1 = spdiags([u u],[1 0],npccx-1,npccx);
               Zx(odd,:) = speye(npccx); %coarse grid points
               Zx(even,:) = Z1;
                             
               Z = kron(Zx,Zx); %2D Operator

        end
end
        
end

