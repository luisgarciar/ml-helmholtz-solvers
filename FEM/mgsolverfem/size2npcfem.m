function [npc,npf] = size2npcfem(s,dim)
%% SIZE2NPCFEM Finds the number of points on a fine coarse grid after standard
%         coarsening given the size of a fine matrix (1D or 2D)
% Input: 
%        s: matrix size
%        d: dimension
%        
% Output:
%        npf, npc: number of points on fine and coarse grids
%
%% 
assert(mod(s,2)==0,'size must be even');

switch dim
    case 1
        npf = s;
    case 2
        npf = sqrt(s);
        assert(floor(npf) == npf,'grid of dim 2 must be square');
end

npc = npf/2;

end

