function [npf,lev] = fem_npc_to_npf(npc,k,par)
%% NPC_TO_NPF Computes the number of points on the finest grid given the
%         number of points on the coarsest grid and discretization
%         requirements
%
% Use:  [npf,lev] = npc_to_npf(npcc,k,ppw)
%
% Input
%       npcc : number of points on coarsest grid (only interior points)
%       k    : wavenumber (for Helmholtz eqn)
%       ppw  : minimum number of points per wavelength   
%
% Output
%       npf : number of interior points on the finest grid
%       l   : number of levels (grids) from coarsest to finest
%
%  Author:      Luis Garcia Ramos, 
%               Institut fur Mathematik, TU Berlin
%               Version 1.0, Jun 2016

%%  
%Case 1: The number of fine grid points grows linearly with the wavenumber
%In this case par equals the number of gridpoints per wavelength
%The number of grid levels lev and points of the fine grid npf is
%chosen according to the rule (2*pi/k*npf) approx par=ppw
if par >1
    ppw = par;
    m   = ppw*k/(2*pi*(npc+1));
    lev = ceil(log2(m));
    lev = max(lev,1);   %at least 1 level
    npf = npc*(2^(lev-1));
    npf = npf-1;
else 
%Case 2: When  par <1, the number of fine grid points grows 
%linearly with respect to k^3/2. , 
%and the gridsize h = 1/(npf+1) satisfies npf approx ceil(k^(3/2))

    m  = (sqrt(k)^3)*(1/(npc+1));
   lev = ceil(log2(m))+1;
   lev = max(lev,1);   %at least 1 level
   npf = (npc+1)*(2^(lev-1))-1;
end

end

