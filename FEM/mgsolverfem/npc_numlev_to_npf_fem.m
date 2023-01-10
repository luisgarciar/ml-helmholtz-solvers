function [npf]= npc_numlev_to_npf_fem(npc,numlev)
%% NPC_NUMLEV_TO_NPF Computes the number of points in a 1D fine grid
%  given the number of points in the coarsest grid and the number 
%  of levels (intermediate grids)
% 
%  Use:
%  [npf] = npc_numlev_to_npf_fem(npc,numlev)
%
%  Input: 
%  npc:    number of points on coarsest grid (1D)
%          (x=0 and interior points)
%
%  numlev: number of levels (grids) including coarsest and finest
%
%  Output:
%  npf:   number of points on finest grid (1D)
%
% Important: 
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 2.0, Oct 2016

npf = 2^(numlev-1)*(npc);

end

