function [npc,npf] = size2npc(s,dim)
%SIZE2NPC Finds the number of points on a fine coarse grid after standard
%         coarsening given the size of a fine matrix (1D or 2D)
%
% Input: 
%        s: matrix size
%        d: dimension
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
end

npc = (npf-1)/2;

end

