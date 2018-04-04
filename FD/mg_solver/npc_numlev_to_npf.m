function [npf]= npc_numlev_to_npf(npc,numlev)
%% NPC_NUMLEV_TO_NPF Computes the number of points in a fine grid
% given the number of points in the coarsest grid and number of levels
% 
%  Use:
% [npf] = npc_numlev_to_npf(npc,numlev)
%
%  Input: 
%  npc:    number of points on coarsest grid (1D)
%
%  numlev: number of levels (grids) including coarsest and finest
%
%  Output:
%  npf:   number of points on finest grid (1D)
%         (x=0 and interior points) 
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 2.0, Oct 2016

%%
%assert(npc > 0,'number of coarse grid points must be positive')
%assert(numlev >0,'number of levels must be positive')

npf = 2^(numlev-1)*(npc+1)-1;

end

