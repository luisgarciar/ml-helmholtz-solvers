function [npc,npf] = size2npc(s,dim,bc)
%SIZE2NPC Finds the number of points on a coarse grid after standard
%             coarsening given the size of a fine matrix (1D or 2D)
%
% Input: 
%        s:  matrix size
%        d:  dimension
%        bc: boundary conditions
%        
% Output:
%        npf, npc: number of points on fine and coarse grids
%
%% 
assert(mod(s,2)==1,'size must be odd');

switch dim
    case 1
        npf = s;
    case 2
         npf = sqrt(s);
         assert(floor(npf) == npf,'grid of dim 2 must be square');

         if strcmp(bc,'som')==1
            npf = npf-2; 
         end
end

npc = (npf-1)/2;

end

